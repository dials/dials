from __future__ import annotations

import copy
import json
import logging
import sys
from collections import OrderedDict

import numpy as np
import scipy.spatial.distance as ssd
from scipy.cluster import hierarchy
from sklearn.cluster import OPTICS

import iotbx.phil
from dxtbx.model import ExperimentList
from libtbx import Auto
from libtbx.phil import scope_extract
from scitbx.array_family import flex

from dials.algorithms.correlation.cluster import ClusterInfo
from dials.algorithms.correlation.plots import (
    linkage_matrix_to_dict,
    plot_dims,
    plot_pca_coords,
    plot_reachability,
    to_plotly_json,
)
from dials.algorithms.symmetry.cosym import CosymAnalysis, cosym_scope
from dials.algorithms.symmetry.cosym.plots import plot_coords, plot_rij_histogram
from dials.array_family.flex import reflection_table
from dials.util import tabulate
from dials.util.exclude_images import get_selection_for_valid_image_ranges
from dials.util.filter_reflections import filtered_arrays_from_experiments_reflections
from dials.util.multi_dataset_handling import select_datasets_on_identifiers

logger = logging.getLogger("dials.algorithms.correlation.analysis")

phil_scope = iotbx.phil.parse(
    """\
partiality_threshold = 0.4
  .type = float(value_min=0, value_max=1)
  .help = "Use reflections with a partiality greater than the threshold."

%s

dimensionality_assessment {
  outlier_rejection = True
    .type = bool
    .help = "Use outlier rejection when determining optimal dimensions for analysis."
  maximum_dimensions = 50
    .type = int
    .help = "Maximum number of dimensions to test for reasonable processing time"
}

significant_clusters {
  min_points_buffer = 0.5
    .type = float(value_min=0, value_max=1)
    .help = "Buffer for minimum number of points required for a cluster in OPTICS algorithm: min_points=(number_of_datasets/number_of_dimensions)*buffer"
  xi = 0.05
    .type = float(value_min=0, value_max=1)
    .help = "xi parameter to determine min steepness to define cluster boundary"
}

hierarchical_clustering {
  linkage_method = *ward average
    .type = choice
    .help = "Linkage method for constructing dendorgrams for both correlation and cosine angle clustering methods"
}
"""
    % cosym_scope,
    process_includes=True,
)
phil_overrides = phil_scope.fetch(
    source=iotbx.phil.parse(
        """\
cc_weights=sigma
weights=count
"""
    )
)
working_phil = phil_scope.fetch(sources=[phil_overrides])


class CorrelationMatrix:
    def __init__(
        self,
        experiments: ExperimentList,
        reflections: list[reflection_table],
        params: scope_extract = None,
        ids_to_identifiers_map: dict = None,
    ):
        """
        Set up the required cosym preparations for determining the correlation matricies
        of a series of input experiments and reflections.

        Args:
            experiments (dxtbx_model_ext.ExperimentList): list of experiments in dials format
            reflections (list): list of dials_array_family_flex_ext.reflection_table objects associated with the experiments
            params (libtbx.phil.scope_extract):
        """

        # Set up experiments, reflections and params

        if params is None:
            params = phil_scope.extract()
        self.params = params

        # if all datasets have been through scaling, a decision about error models has
        # been made, so don't apply any further sigma correction
        apply_sigma_correction = not all(s for s in experiments.scaling_models())

        self._reflections = []
        self.ids_to_identifiers_map = ids_to_identifiers_map

        if len(reflections) == len(experiments):
            for refl, expt in zip(reflections, experiments):
                sel = get_selection_for_valid_image_ranges(refl, expt)
                self._reflections.append(refl.select(sel))
        else:
            sys.exit(
                "Number of reflection tables does not match number of experiments."
            )

        # Initial filtering to remove experiments and reflections that do not meet the minimum number of reflections required (params.min_reflections)

        self._experiments, self._reflections = self._filter_min_reflections(
            experiments, self._reflections
        )

        # Used for optional json creation that is in a format friendly for import and analysis (for future development)
        # Also to retain dataset ids when used in multiplex
        if self.ids_to_identifiers_map is None:
            self.ids_to_identifiers_map = {}
            for table in self._reflections:
                self.ids_to_identifiers_map.update(table.experiment_identifiers())

        self.labels = list(dict.fromkeys(self.ids_to_identifiers_map))
        self._labels_all = flex.size_t(self.labels)

        # Filter reflections that do not meet partiality threshold or default I/Sig(I) criteria

        datasets = filtered_arrays_from_experiments_reflections(
            self._experiments,
            self._reflections,
            outlier_rejection_after_filter=False,
            partiality_threshold=params.partiality_threshold,
        )

        self.unmerged_datasets = datasets

        # Merge intensities to prepare for cosym analysis

        self.datasets = self._merge_intensities(datasets)

        # Set required params for cosym to skip symmetry determination and reduce dimensions

        self.params.__inject__("lattice_group", self.datasets[0].space_group_info())
        self.params.__inject__("space_group", self.datasets[0].space_group_info())
        self.params.__inject__("lattice_symmetry_max_delta", 0.0)
        self.params.__inject__("best_monoclinic_beta", True)

        self.cosym_analysis = CosymAnalysis(
            self.datasets, self.params, apply_sigma_correction=apply_sigma_correction
        )

    def _merge_intensities(self, datasets: list) -> list:
        """
        Merge intensities and elimate systematically absent reflections.

        Args:
            datasets(list): list of cctbx.miller.array objects
        Returns:
            datasets_sys_absent_eliminated(list): list of merged cctbx.miller.array objects
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

    def _filter_min_reflections(
        self, experiments: ExperimentList, reflections: list[reflection_table]
    ) -> tuple[ExperimentList, list[reflection_table]]:
        """
        Filter all datasets that have less than the specified number of reflections.

        Args:
            experiments (dxtbx_model_ext.ExperimentList): list of experiments in dials format
            reflections (list): list of dials_array_family_flex_ext.reflection_table objects associated with the experiments

        Returns:
            filtered_datasets (tuple): tuple of filtered datasets

        """
        identifiers = []

        for expt, refl in zip(experiments, reflections):
            if len(refl) >= self.params.min_reflections:
                identifiers.append(expt.identifier)

        filtered_datasets = select_datasets_on_identifiers(
            experiments, reflections, use_datasets=identifiers
        )

        return filtered_datasets

    def calculate_matrices(self):
        """
        Runs the required algorithms within dials.cosym to calculate the rij matrix and optimise the coordinates.
        These results are passed into matrix computation functions to calculate the cosine similarity (cos-angle) and correlation matrices
        as well as the corresponding clustering.
        """

        # Cosym algorithm to calculate the rij matrix (CC matrix when symmetry known)
        self.cosym_analysis._intialise_target()

        # Cosym proceedures to calculate the cos-angle matrix
        if (
            len(self.unmerged_datasets)
            <= self.params.dimensionality_assessment.maximum_dimensions
        ):
            dims_to_test = len(self.unmerged_datasets)
        else:
            dims_to_test = self.params.dimensionality_assessment.maximum_dimensions

        self._dimension_optimisation_data = {}

        if self.params.dimensions is Auto:
            (
                self._dimension_optimisation_data["dimensions"],
                self._dimension_optimisation_data["functional"],
            ) = self.cosym_analysis._determine_dimensions(
                dims_to_test,
                outlier_rejection=self.params.dimensionality_assessment.outlier_rejection,
            )
        self.cosym_analysis._optimise(
            self.cosym_analysis.params.minimization.engine,
            max_iterations=self.cosym_analysis.params.minimization.max_iterations,
            max_calls=self.cosym_analysis.params.minimization.max_calls,
        )

        self.cosym_analysis._principal_component_analysis(cluster=True)

        # Convert cosym output into cc/cos matrices in correct form and compute linkage matrices as clustering method
        (
            self.correlation_matrix,
            self.cc_linkage_matrix,
        ) = self.compute_correlation_coefficient_matrix(
            self.cosym_analysis.target.rij_matrix,
            self.params.hierarchical_clustering.linkage_method,
        )

        self.correlation_clusters = self.cluster_info(
            linkage_matrix_to_dict(self.cc_linkage_matrix)
        )

        logger.info("\nIntensity correlation clustering summary:")
        self.cc_table = ClusterInfo.as_table(self.correlation_clusters)
        logger.info(tabulate(self.cc_table, headers="firstrow", tablefmt="rst"))
        self.cos_angle, self.cos_linkage_matrix = self.compute_cos_angle_matrix(
            self.cosym_analysis.coords,
            self.params.hierarchical_clustering.linkage_method,
        )

        self.cos_angle_clusters = self.cluster_info(
            linkage_matrix_to_dict(self.cos_linkage_matrix)
        )

        logger.info("\nCos(angle) clustering summary:")
        self.cos_table = ClusterInfo.as_table(self.cos_angle_clusters)
        logger.info(tabulate(self.cos_table, headers="firstrow", tablefmt="rst"))

        logger.info("\nEvaluating Significant Clusters from Cosine-Angle Coordinates:")
        self.cluster_cosine_coords()

    @staticmethod
    def compute_correlation_coefficient_matrix(
        correlation_matrix: np.ndarray,
        linkage_method: str,
    ) -> tuple(np.ndarray, np.ndarray):
        """
        Computes the correlation matrix and clustering linkage matrix from the rij cosym matrix.

        Args:
            correlation_matrix(numpy.ndarray): pair-wise matrix of correlation coefficients
            linkage_method(str): linkage method for hierarchical clustering

        Returns:
            correlation_matrix(numpy.ndarray): correlation matrix with corrections to diagonals and accounting for floating point errors
            cc_linkage_matrix(numpy.ndarray): linkage matrix describing clustering of correlation matrix in dendrogram-style

        """

        logger.info("\nCalculating Correlation Matrix (rij matrix - see dials.cosym)")
        logger.info(
            f"Hierarchical clustering performed using the {linkage_method} linkage method"
        )

        # Make diagonals equal to 1 (each dataset correlated with itself)
        np.fill_diagonal(correlation_matrix, 1)

        # clip values of correlation matrix to account for floating point errors
        correlation_matrix.clip(-1, 1, out=correlation_matrix)

        # Convert to distance matrix rather than correlation
        diffraction_dissimilarity = 1 - correlation_matrix
        try:
            assert ssd.is_valid_dm(diffraction_dissimilarity, tol=1e-12)
        except AssertionError:
            sys.exit(
                "Correlation matrix does not give a valid distance matrix. Distance matrix is either non-symmetric or does not have a zero-diagonal."
            )

        # convert the redundant n*n square matrix form into a condensed nC2 array
        cc_dist_mat = ssd.squareform(diffraction_dissimilarity, checks=False)

        # Clustering method
        cc_linkage_matrix = hierarchy.linkage(cc_dist_mat, method=linkage_method)

        return correlation_matrix, cc_linkage_matrix

    @staticmethod
    def compute_cos_angle_matrix(
        coords: np.ndarray,
        linkage_method: str,
    ) -> tuple(np.ndarray, np.ndarray):
        """
        Computes the cos_angle matrix and clustering linkage matrix from the optimized cosym coordinates.
        Args:
            coords(numpy.ndarray): matrix of coordinates output from cosym optimisation
            linkage_method(str): linkage method for hierarchical clustering

        Returns:
            cos_angle(numpy.ndarray): pair-wise cos angle matrix
            cos_linkage_matrix(numpy.ndarray): linkage matrix describing clustering of cos angle matrix in dendrogram-style

        """
        logger.info(
            "\nCalculating Cos Angle Matrix from optimised cosym coordinates (see dials.cosym)"
        )
        logger.info(
            f"Hierarchical clustering performed using the {linkage_method} linkage method"
        )

        # Convert coordinates to cosine distances and then reversed so closer cosine distances have higher values to match CC matrix
        cos_dist_mat = ssd.pdist(coords, metric="cosine")
        cos_angle = 1 - ssd.squareform(cos_dist_mat)

        # Clustering method
        cos_linkage_matrix = hierarchy.linkage(cos_dist_mat, method=linkage_method)

        return cos_angle, cos_linkage_matrix

    def cluster_cosine_coords(self):
        """
        Cluster cosine coords using the OPTICS algorithm to determine significant clusters.
        """

        logger.info("Using OPTICS Algorithm (M. Ankerst et al, 1999, ACM SIGMOD)")

        # Minimum number required to make sense for such algorithms

        min_points = max(
            5,
            int(
                (len(self.unmerged_datasets) / self.cosym_analysis.target.dim)
                * self.params.significant_clusters.min_points_buffer
            ),
        )

        # Check for very small datasets
        if len(self.unmerged_datasets) < min_points:
            min_points = len(self.unmerged_datasets)
            logger.info(
                "WARNING: less than 5 samples present, OPTICS not optimised for very small datasets."
            )

        logger.info(f"Setting Minimum Samples to {min_points}")

        # Fit OPTICS model and determine number of clusters

        optics_model = OPTICS(
            min_samples=min_points, xi=self.params.significant_clusters.xi
        )

        optics_model.fit(self.cosym_analysis.coords)

        self.cluster_labels = optics_model.labels_

        # Reachability plot data

        self.optics_reachability = optics_model.reachability_[optics_model.ordering_]
        self.optics_reachability_labels = optics_model.labels_[optics_model.ordering_]

        # Match each dataset to an OPTICS cluster and make them Cluster Objects

        unique_labels = sorted(dict.fromkeys(self.cluster_labels))

        if -1 in unique_labels:
            unique_labels.remove(-1)

        outliers = 0
        for i in self.cluster_labels:
            if i == -1:
                outliers += 1

        logger.info(
            f"OPTICS identified {len(unique_labels)} clusters and {outliers} outlier datasets."
        )

        sig_cluster_dict = OrderedDict()

        # Seems complex to deal with cases where datasets have been filtered
        # and fit into how other code already set up

        for_cluster_ids = list(range(1, len(self.labels) + 1))

        # Set up number of clusters
        for i in unique_labels:
            ids = []
            # Iterate through cluster labels from OPTICS to match to IDS
            for main_idx, clust in enumerate(self.cluster_labels):
                if i == clust:
                    # Pre-existing code requieres +1 to everything
                    for idx1, val1 in enumerate(for_cluster_ids):
                        # Real labels start at 0
                        for idx2, val2 in enumerate(self.labels):
                            # Need to match index due to some datasets possibly being filtered out
                            if idx1 == idx2 and idx1 == main_idx:
                                ids.append(val1)
            cluster_dict = {"datasets": ids}
            sig_cluster_dict[i] = cluster_dict

        self.significant_clusters = self.cluster_info(sig_cluster_dict)

        for i in self.significant_clusters:
            logger.info(i)

    def cluster_info(self, cluster_dict: dict) -> list:
        """
        Generate list of cluster objects with associated statistics.
        Args:
            cluster_dict(dict): dictionary of clusters (generated from linkage_matrix_to_dict)

        Returns:
            info(list): list of ClusterInfo objects to describe all clusters of a certain type (ie correlation or cos angle)
        """

        info = []
        for cluster_id, cluster in cluster_dict.items():
            uc_params = [flex.double() for i in range(6)]
            for j in cluster["datasets"]:
                uc_j = self.datasets[j - 1].unit_cell().parameters()
                for i in range(6):
                    uc_params[i].append(uc_j[i])
            average_uc = [flex.mean(uc_params[i]) for i in range(6)]
            intensities_cluster = []
            labels_cluster = []
            ids = [self._labels_all[id - 1] for id in cluster["datasets"]]
            for idx, k in zip(self._labels_all, self.unmerged_datasets):
                if idx in ids:
                    intensities_cluster.append(k)
                    labels_cluster.append(idx)
            merged = None
            for d in intensities_cluster:
                if merged is None:
                    merged = copy.deepcopy(d)
                else:
                    merged = merged.concatenate(d, assert_is_similar_symmetry=False)
            merging = merged.merge_equivalents()
            merged_intensities = merging.array()
            multiplicities = merging.redundancies()
            info.append(
                ClusterInfo(
                    cluster_id,
                    labels_cluster,
                    flex.mean(multiplicities.data().as_double()),
                    merged_intensities.completeness(),
                    unit_cell=average_uc,
                    height=cluster.get("height"),
                )
            )
        return info

    def convert_to_html_json(self):
        """
        Prepares the required dataset tables and converts analysis into the required format for HTML output.
        """

        # Convert the cosine and cc matrices into a plotly json format for output graphs

        self.cc_json = to_plotly_json(
            self.correlation_matrix,
            self.cc_linkage_matrix,
            labels=self.labels,
            matrix_type="correlation",
        )
        self.cos_json = to_plotly_json(
            self.cos_angle,
            self.cos_linkage_matrix,
            labels=self.labels,
            matrix_type="cos_angle",
        )

        self.rij_graphs = OrderedDict()

        self.rij_graphs.update(
            plot_rij_histogram(self.correlation_matrix, key="cosym_rij_histogram_sg")
        )

        if self._dimension_optimisation_data:
            self.rij_graphs.update(
                plot_dims(
                    self._dimension_optimisation_data["dimensions"],
                    self._dimension_optimisation_data["functional"],
                )
            )

        self.rij_graphs.update(
            plot_reachability(
                np.arange(len(self.optics_reachability)),
                self.optics_reachability,
                self.optics_reachability_labels,
            )
        )

        cluster_labels = self.cluster_labels
        axes_labels = {
            str(i): f"PC {i + 1} ({var:.1f}%)"
            for i, var in enumerate(self.cosym_analysis.explained_variance_ratio * 100)
        }

        # Cosym coordinates are aligned with the principal components
        # Note that the reduced coordinates are NOT used because want to retain information
        # about length/angles from the original cosym procedure
        # Instead, simply rotate the data to align with the eigenvectors

        cob_mat = np.linalg.inv(np.transpose(self.cosym_analysis.pca_components))

        coords = np.array([])
        for data in self.cosym_analysis.coords:
            rot_coord = np.matmul(cob_mat, data)
            if coords.size == 0:
                coords = rot_coord
            else:
                coords = np.vstack([coords, rot_coord])

        if coords.shape[1] > 6:
            coords = coords[:, 0:6]
            dim_list = list(range(0, 6))
        else:
            dim_list = list(range(0, self.cosym_analysis.target.dim))

        self.pca_plot = plot_pca_coords(
            coords,
            axes_labels,
            cluster_labels,
            dim_list,
        )

        # Separately also plot the two most significant principal components for ease of visualisation

        self.rij_graphs.update(
            plot_coords(
                coords[:, 0:2],
                self.cluster_labels,
                key="cosym_coordinates_principal_components",
            )
        )

        updated_help = """\
The outcome of the cosym multi-dimensional analysis rotated by the components identified using
Principal Component Analysis and projected in 2D. Each point corresponds to an individual data set.
The lengths of the vectors are inversely related to the amount of random error in each dataset;
however, because the plot is a projection in 2D, the coordinates in this graph are not directly
representative of the length of the vector in multi-dimensional space (unless the analysis
is 2-dimensional). The angular separation between any pair, or groups, of vectors is a measure of the
systematic differences between the data sets (could be the presence of non-isomorphism). This projection
summarises the effect of the first two principal components.
        """

        self.rij_graphs["cosym_coordinates_principal_components"]["help"] = updated_help
        self.rij_graphs["cosym_coordinates_principal_components"]["layout"]["title"] = (
            "Cosym Coordinates Rotated by Principal Components"
        )
        self.rij_graphs["cosym_coordinates_principal_components"]["layout"]["xaxis"][
            "title"
        ] = axes_labels["0"]
        self.rij_graphs["cosym_coordinates_principal_components"]["layout"]["yaxis"][
            "title"
        ] = axes_labels["1"]

        # Generate the table for the html that lists all datasets and image paths present in the analysis

        paths = enumerate(e.imageset.paths()[0] for e in self._experiments)
        self.table_list = [["Experiment/Image Number", "Image Path"], *map(list, paths)]

    def convert_to_importable_json(self, linkage_matrix: np.ndarray) -> OrderedDict:
        """
        Generate a json file of the linkage matrices with unique identifiers rather than dataset numbers
        May be useful for future developments to import this for further clustering analysis.
        Args:
            linkage_matrix(numpy.ndarray): linkage matrix from hierarchy.linkage methods
        Returns:
            linkage_mat_as_dict(collections.OrderedDict): linkage matrix converted to dictionary with datasets replaced with dials unique identifiers
        """
        linkage_mat_as_dict = linkage_matrix_to_dict(linkage_matrix)
        for d in linkage_mat_as_dict.values():
            # Difference in indexing between linkage_mat_as_dict and datasets, so have i-1
            real_num = [self.labels[i - 1] for i in d["datasets"]]
            d["datasets"] = [self.ids_to_identifiers_map[i] for i in real_num]

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

        with open(self.params.output.json, "w") as f:
            json.dump(combined_json_dict, f)

    def output_clusters(self):
        """
        Output dials .expt/.refl files for each significant cluster identified in the OPTICS analysis of the cosine angle clustering.
        """

        from dials.array_family import flex

        for i in self.significant_clusters:
            expts = ExperimentList()
            refls = []
            for idx in i.labels:
                expts.append(self._experiments[idx])
                refls.append(self._reflections[idx])
            joint_refl = flex.reflection_table.concat(refls)
            expts.as_file(f"cluster_{i.cluster_id}.expt")
            joint_refl.as_file(f"cluster_{i.cluster_id}.refl")

        logger.info("Identified clusters output as .expt/.refl files.")
