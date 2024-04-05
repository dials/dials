from __future__ import annotations

import json
import logging
import sys
from collections import OrderedDict

import numpy as np
import scipy.spatial.distance as ssd
from scipy.cluster import hierarchy

import iotbx.phil

from dials.algorithms.correlation.plots import linkage_matrix_to_dict, to_plotly_json
from dials.algorithms.symmetry.cosym import CosymAnalysis
from dials.algorithms.symmetry.cosym.plots import plot_coords, plot_rij_histogram
from dials.util.exclude_images import get_selection_for_valid_image_ranges
from dials.util.filter_reflections import filtered_arrays_from_experiments_reflections
from dials.util.multi_dataset_handling import select_datasets_on_identifiers

logger = logging.getLogger("dials.algorithms.correlation.analysis")

phil_scope = iotbx.phil.parse(
    """\
partiality_threshold = 0.4
  .type = float(value_min=0, value_max=1)
  .help = "Use reflections with a partiality greater than the threshold."

include scope dials.algorithms.symmetry.cosym.phil_scope

relative_length_tolerance = 0.05
  .type = float(value_min=0)
  .help = "Datasets with unit cell lengths are only accepted if within this relative tolerance of the median cell."

absolute_angle_tolerance = 2
  .type = float(value_min=0)
  .help = "Datasets with unit cell angles are only accepted if within this absolute tolerance of the median cell."

min_reflections = 10
  .type = int(value_min=1)
  .help = "The minimum number of reflections per experiment."
""",
    process_includes=True,
)


class CorrelationMatrix:
    def __init__(self, experiments, reflections, params=None):
        """
        Set up the required cosym preparations for determining the correlation matricies
        of a series of input experiments and reflections.
        """

        # Set up experiments, reflections and params

        if params is None:
            params = phil_scope.extract()
        self.params = params
        self._reflections = []

        if len(reflections) == len(experiments):
            for refl, expt in zip(reflections, experiments):
                sel = get_selection_for_valid_image_ranges(refl, expt)
                self._reflections.append(refl.select(sel))
        else:
            sys.exit("Number of reflections and experiments do not match.")

        # Initial filtering to remove experiments and reflections that do not meet the minimum number of reflections required (params.min_reflections)

        self._experiments, self._reflections = self._filter_min_reflections(
            experiments, self._reflections
        )

        # Used for optional json creation that is in a format friendly for import and analysis (for future development)
        self.ids_to_identifiers_map = {}
        for table in self._reflections:
            self.ids_to_identifiers_map.update(table.experiment_identifiers())

        # Filter reflections that do not meet partiality threshold or default I/Sig(I) criteria

        datasets = filtered_arrays_from_experiments_reflections(
            self._experiments,
            self._reflections,
            outlier_rejection_after_filter=False,
            partiality_threshold=params.partiality_threshold,
        )

        # Merge intensities to prepare for cosym analysis

        self.datasets = self._merge_intensities(datasets)

        # Set required params for cosym to skip symmetry determination and reduce dimensions

        self.params.lattice_group = self.datasets[0].space_group_info()
        self.params.space_group = self.datasets[0].space_group_info()

        self.cosym_analysis = CosymAnalysis(self.datasets, self.params)

    def _merge_intensities(self, datasets):
        """
        Merge intensities and elimate systematically absent reflections
        """

        individual_merged_intensities = []
        for unmerged in datasets:
            individual_merged_intensities.append(
                unmerged.merge_equivalents().array().set_info(unmerged.info())
            )
        datasets_sys_absent_eliminated = [
            d.eliminate_sys_absent(integral_only=True).primitive_setting()
            for d in individual_merged_intensities
        ]

        return datasets_sys_absent_eliminated

    def _filter_min_reflections(self, experiments, reflections):
        """
        Filter all datasets that have less than the specified number of reflections.
        """
        identifiers = []

        for expt, refl in zip(experiments, reflections):
            if len(refl) >= self.params.min_reflections:
                identifiers.append(expt.identifier)

        return select_datasets_on_identifiers(
            experiments, reflections, use_datasets=identifiers
        )

    def calculate_matrices(self):
        """
        Runs the required algorithms within dials.cosym to calculate the rij matrix and optimise the coordinates.
        These results are passed into matrix computation functions to calculate the cos-angle and correlation matrices
        and corresponding clustering.
        """

        # Cosym algorithm to calculate the rij matrix (CC matrix when symmetry known)
        self.cosym_analysis._intialise_target()

        # Cosym proceedures to calculate the cos-angle matrix
        self.cosym_analysis._determine_dimensions()
        self.cosym_analysis._optimise(
            self.cosym_analysis.params.minimization.engine,
            max_iterations=self.cosym_analysis.params.minimization.max_iterations,
            max_calls=self.cosym_analysis.params.minimization.max_calls,
        )

        # Convert cosym output into cc/cos matrices in correct form and compute linkage matrices as clustering method
        (
            self.correlation_matrix,
            self.cc_linkage_matrix,
        ) = self.compute_correlation_coefficient_matrix(
            self.cosym_analysis.target.rij_matrix
        )
        self.cos_angle, self.cos_linkage_matrix = self.compute_cos_angle_matrix(
            self.cosym_analysis.coords
        )

    @staticmethod
    def compute_correlation_coefficient_matrix(correlation_matrix):
        """
        Computes the correlation matrix and clustering linkage matrix from the rij cosym matrix.
        """

        logger.info("\nCalculating Correlation Matrix (rij matrix - see dials.cosym)\n")

        # Make diagonals equal to 1 (each dataset correlated with itself)
        for i in range(correlation_matrix.shape[0]):
            correlation_matrix[i, i] = 1

        # clip values of correlation matrix to account for floating point errors
        correlation_matrix[np.where(correlation_matrix < -1)] = -1
        correlation_matrix[np.where(correlation_matrix > 1)] = 1

        # Convert to distance matrix rather than correlation
        diffraction_dissimilarity = 1 - correlation_matrix
        assert ssd.is_valid_dm(diffraction_dissimilarity, tol=1e-12)

        # convert the redundant n*n square matrix form into a condensed nC2 array
        cc_dist_mat = ssd.squareform(diffraction_dissimilarity, checks=False)

        # Clustering method
        cc_linkage_matrix = hierarchy.linkage(cc_dist_mat, method="ward")

        return correlation_matrix, cc_linkage_matrix

    @staticmethod
    def compute_cos_angle_matrix(coords):
        """
        Computes the cos_angle matrix and clustering linkage matrix from the optimized cosym coordinates.
        """

        logger.info(
            "Calculating Cos Angle Matrix from optimised cosym coordinates (see dials.cosym)\n"
        )

        # Convert coordinates to cosine distances and then reversed so closer cosine distances have higher values to match CC matrix
        cos_dist_mat = ssd.pdist(coords, metric="cosine")
        cos_angle = 1 - ssd.squareform(cos_dist_mat)

        # Clustering method
        cos_linkage_matrix = hierarchy.linkage(cos_dist_mat, method="ward")

        return cos_angle, cos_linkage_matrix

    def convert_to_html_json(self):
        """
        Prepares the required dataset tables and converts analysis into the required format for HTML output.
        """

        # Convert the cosine and cc matrices into a plotly json format for output graphs

        labels = list(range(0, len(self._experiments)))

        self.cc_json = to_plotly_json(
            self.correlation_matrix,
            self.cc_linkage_matrix,
            labels=labels,
            matrix_type="correlation",
        )
        self.cos_json = to_plotly_json(
            self.cos_angle,
            self.cos_linkage_matrix,
            labels=labels,
            matrix_type="cos_angle",
        )

        self.rij_graphs = OrderedDict()

        self.rij_graphs.update(
            plot_rij_histogram(self.correlation_matrix, key="cosym_rij_histogram_sg")
        )

        self.rij_graphs.update(
            plot_coords(self.cosym_analysis.coords, key="cosym_coordinates_sg")
        )

        # Generate the table for the html that lists all datasets and image paths present in the analysis

        path_list = []
        self.table_list = [["Experiment/Image Number", "Image Path"]]

        for i in self._experiments:
            path_list.append(i.imageset.paths()[0])

        ids = list(range(0, len(path_list)))

        for i, j in zip(ids, path_list):
            self.table_list.append([i, j])

    def convert_to_importable_json(self, linkage_matrix):
        """
        Generate a json file of the linkage matrices with unique identifiers rather than dataset numbers
        May be useful for future developments to import this for further clustering analysis
        """

        linkage_mat_as_dict = linkage_matrix_to_dict(linkage_matrix)
        for i in linkage_mat_as_dict:
            # Difference in indexing between linkage_mat_as_dict and datasets, so have i-1
            old_datasets = [i - 1 for i in linkage_mat_as_dict[i]["datasets"]]
            datasets_as_ids = [self.ids_to_identifiers_map[j] for j in old_datasets]
            linkage_mat_as_dict[i]["datasets"] = datasets_as_ids

        return linkage_mat_as_dict

    def output_json(self):
        """
        Outputs the cos and cc json files containing correlation/cos matrices and linkage details.
        """

        linkage_mat_as_dict_cc = self.convert_to_importable_json(self.cc_linkage_matrix)
        linkage_mat_as_dict_cos = self.convert_to_importable_json(
            self.cos_linkage_matrix
        )

        combined_json_dict = {
            "correlation_matrix_clustering": linkage_mat_as_dict_cc,
            "correlation_matrix": self.correlation_matrix.tolist(),
            "cos_matrix_clustering": linkage_mat_as_dict_cos,
            "cos_angle_matrix": self.cos_angle.tolist(),
        }

        json_str = json.dumps(combined_json_dict)
        with open(self.params.output.json, "w") as f:
            f.write(json_str)
