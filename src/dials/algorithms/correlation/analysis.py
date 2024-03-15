from __future__ import annotations

import json
import logging
from collections import OrderedDict

import numpy as np
import scipy.spatial.distance as ssd
from scipy.cluster import hierarchy

import iotbx.phil

from dials.algorithms.correlation.plots import to_plotly_json
from dials.algorithms.symmetry.cosym import CosymAnalysis
from dials.algorithms.symmetry.cosym.plots import plot_coords, plot_rij_histogram
from dials.util.exclude_images import get_selection_for_valid_image_ranges
from dials.util.filter_reflections import filtered_arrays_from_experiments_reflections
from dials.util.multi_dataset_handling import select_datasets_on_identifiers
from dials.util.observer import Subject

logger = logging.getLogger("dials.algorithms.correlation.analysis")

phil_scope = iotbx.phil.parse(
    """\
partiality_threshold = 0.4
  .type = float
  .help = "Use reflections with a partiality above the threshold."

include scope dials.algorithms.symmetry.cosym.phil_scope

relative_length_tolerance = 0.05
  .type = float(value_min=0)

absolute_angle_tolerance = 2
  .type = float(value_min=0)

min_reflections = 10
  .type = int(value_min=1)
  .help = "The minimum number of reflections per experiment."
""",
    process_includes=True,
)


class CorrelationMatrix(Subject):
    def __init__(self, experiments, reflections, params=None):
        """
        Set up the required cosym preparations for determining the correlation matricies
        of a series of input experiments and reflections.
        Args:
          experiments (dxtbx_model_ext.ExperimentList): dials experiments
          reflections (list): dials reflections
          params (libtbx.phil.scope_extract): experimental parameters
        """

        if params is None:
            params = phil_scope.extract()
        self.params = params

        self._reflections = []
        for refl, expt in zip(reflections, experiments):
            sel = get_selection_for_valid_image_ranges(refl, expt)
            self._reflections.append(refl.select(sel))

        self._experiments, self._reflections = self._filter_min_reflections(
            experiments, self._reflections
        )
        self.ids_to_identifiers_map = {}
        for table in self._reflections:
            self.ids_to_identifiers_map.update(table.experiment_identifiers())
        self.identifiers_to_ids_map = {
            value: key for key, value in self.ids_to_identifiers_map.items()
        }

        if "inverse_scale_factor" in reflections[0]:
            # FIXME - can we also just use the filtered_arrays_from_experiments_reflections function?
            from dials.report.analysis import scaled_data_as_miller_array

            datasets = []
            for expt, r in zip(self._experiments, self._reflections):
                sel = ~r.get_flags(r.flags.bad_for_scaling, all=False)
                sel &= r["inverse_scale_factor"] > 0
                datasets.append(scaled_data_as_miller_array([r], [expt]))

        else:
            datasets = filtered_arrays_from_experiments_reflections(
                self._experiments,
                self._reflections,
                outlier_rejection_after_filter=False,
                partiality_threshold=params.partiality_threshold,
            )
        individual_merged_intensities = []
        for unmerged in datasets:
            individual_merged_intensities.append(
                unmerged.merge_equivalents().array().set_info(unmerged.info())
            )
        self.datasets = [
            d.eliminate_sys_absent(integral_only=True).primitive_setting()
            for d in individual_merged_intensities
        ]

        self.params.lattice_group = self.datasets[0].space_group_info()
        self.params.space_group = self.datasets[0].space_group_info()

        self.cosym_analysis = CosymAnalysis(self.datasets, self.params)

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
        Calculates the cc and cos angle correlation matrices for the input dataset.
        The cc matrix is converted to a distance matrix, and both are condensed into a nC2 array.
        Matrices are converted to plotly json files for visual output.
        """

        self.cosym_analysis._intialise_target()
        self.cosym_analysis._determine_dimensions()
        self.cosym_analysis._optimise(
            self.cosym_analysis.params.minimization.engine,
            max_iterations=self.cosym_analysis.params.minimization.max_iterations,
            max_calls=self.cosym_analysis.params.minimization.max_calls,
        )

        logger.info("\nCalculating Correlation Matrix (rij matrix - see dials.cosym)\n")

        correlation_matrix = self.cosym_analysis.target.rij_matrix

        for i in range(correlation_matrix.shape[0]):
            correlation_matrix[i, i] = 1

        # clip values of correlation matrix to account for floating point errors
        correlation_matrix[np.where(correlation_matrix < -1)] = -1
        correlation_matrix[np.where(correlation_matrix > 1)] = 1
        diffraction_dissimilarity = 1 - correlation_matrix

        assert ssd.is_valid_dm(diffraction_dissimilarity, tol=1e-12)
        # convert the redundant n*n square matrix form into a condensed nC2 array
        cc_dist_mat = ssd.squareform(diffraction_dissimilarity, checks=False)

        cc_linkage_matrix = hierarchy.linkage(cc_dist_mat, method="average")

        logger.info(
            "Calculating Cos Angle Matrix from optimised cosym coordinates (see dials.cosym)\n"
        )

        cos_dist_mat = ssd.pdist(self.cosym_analysis.coords, metric="cosine")
        cos_angle = 1 - ssd.squareform(cos_dist_mat)
        cos_linkage_matrix = hierarchy.linkage(cos_dist_mat, method="average")

        self.cc_linkage_matrix = cc_linkage_matrix
        self.correlation_matrix = correlation_matrix
        self.cos_angle = cos_angle
        self.cos_linkage_matrix = cos_linkage_matrix

        labels = list(range(0, len(self._experiments)))

        self.cc_json = to_plotly_json(
            correlation_matrix,
            cc_linkage_matrix,
            labels=labels,
            matrix_type="correlation",
        )
        self.cos_json = to_plotly_json(
            cos_angle, cos_linkage_matrix, labels=labels, matrix_type="cos_angle"
        )

        self.rij_graphs = OrderedDict()

        self.rij_graphs.update(
            plot_rij_histogram(correlation_matrix, key="cosym_rij_histogram_sg")
        )

        self.rij_graphs.update(
            plot_coords(self.cosym_analysis.coords, key="cosym_coordinates_sg")
        )

        path_list = []
        self.table_list = [["Experiment/Image Number", "Image Path"]]

        for i in self._experiments:
            j = i.imageset
            path_list.append(j.paths()[0])

        ids = list(range(0, len(path_list)))

        for i, j in zip(ids, path_list):
            self.table_list.append([i, j])

    def output_json(self):
        """
        Outputs the cos and cc json files used for graphing.
        """

        json_str = json.dumps(self.cc_json)
        with open(self.params.output.cc_json, "w") as f:
            f.write(json_str)
        json_str = json.dumps(self.cos_json)
        with open(self.params.output.cos_json, "w") as f:
            f.write(json_str)
