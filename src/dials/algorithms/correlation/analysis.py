from __future__ import annotations

import copy
import itertools
import json
import logging
import math
from collections import OrderedDict

import numpy as np
import scipy.spatial.distance as ssd
from scipy.cluster import hierarchy

import iotbx.phil
from cctbx.miller import binned_data

from dials.algorithms.correlation.plots import to_plotly_json
from dials.algorithms.symmetry.cosym import CosymAnalysis
from dials.algorithms.symmetry.cosym.plots import plot_coords, plot_rij_histogram
from dials.array_family import flex
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

        # self.params.weights = "standard_error"

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

        # DO WE WANT TO ASSERT THIS?!?!?!?!?!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

        logger.info("Applying matrix corrections via scary maths\n")

        # TEMP - MAKE THIS NICER LATER
        self.cc_linkage_matrix = cc_linkage_matrix
        self.correlation_matrix = correlation_matrix
        self.cos_angle = cos_angle
        self.cos_linkage_matrix = cos_linkage_matrix

        cc_weighted_matrix = np.zeros(correlation_matrix.shape)
        for row, i in enumerate(self.datasets):
            for column, j in enumerate(self.datasets):
                m1 = copy.deepcopy(i)
                m2 = copy.deepcopy(j)
                m1_int, m2_int = m1.common_sets(m2)
                wcc, neff = self.weighted_cchalf(m1_int, m2_int)
                cc_weighted_matrix[row, column] = wcc

        difference_matrix = abs(correlation_matrix - cc_weighted_matrix)
        logger.info("Difference between cc matrix and cc weighted matrix")
        logger.info(difference_matrix)
        logger.info("average of difference matrix")
        logger.info(np.average(difference_matrix))

        # DO WE WANT TO ASSERT THIS?!?!?!?!?!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for i in range(cc_weighted_matrix.shape[0]):
            cc_weighted_matrix[i, i] = 1

        # clip values of correlation matrix to account for floating point errors
        cc_weighted_matrix[np.where(cc_weighted_matrix < -1)] = -1
        cc_weighted_matrix[np.where(cc_weighted_matrix > 1)] = 1
        diffraction_dissimilarity_corr = 1 - cc_weighted_matrix

        assert ssd.is_valid_dm(diffraction_dissimilarity_corr, tol=1e-12)
        # convert the redundant n*n square matrix form into a condensed nC2 array
        cc_dist_mat_corr = ssd.squareform(diffraction_dissimilarity_corr, checks=False)

        cc_linkage_matrix_corr = hierarchy.linkage(cc_dist_mat_corr, method="average")

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

        self.cc_corr_json = to_plotly_json(
            cc_weighted_matrix,
            cc_linkage_matrix_corr,
            labels=labels,
            matrix_type="correlation",
        )
        self.cos_corr_json = self.cos_json

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
        json_str = json.dumps(self.cc_corr_json)
        with open(self.params.output.cc_corr_json, "w") as f:
            f.write(json_str)
        json_str = json.dumps(self.cos_corr_json)
        with open(self.params.output.cos_corr_json, "w") as f:
            f.write(json_str)

    def CompleteCovarianceMatrix(self, Corr, S):
        # Corr is the initial (possibly invalid) covariance matrix
        # S is a matrix with 1 or 0 for specified covariances

        C = copy.deepcopy(Corr)

        # number of rows
        n = np.shape(C)[0]
        p = 0.98
        # OCTAVE starts from index 1, python starts from index 0
        for i in range(1, n):
            # Diagonals defined to be known, so this finds all that are not diagonals (and assume triangular symmetry?)
            SS = S[i, 0:i]
            # Makes a np.array with the INDICES that are non-zero - need [0] because gets stored as a weird tuple
            # MAYBE RETHINK THIS LATER
            ind = (SS == 1).nonzero()[0]
            # Makes a np.array with the INDICES that are zero - need [0] because gets stored as a weird tuple
            # MAYBE RETHINK THIS LATER
            nind = (SS == 0).nonzero()[0]
            # Convert specified covariances into a matrix that looks like this:
            # H = [A   B'  x]
            #     [B   D   c]
            #     [x'  c'  d]
            d = C[i, i]
            c = C[ind, i]
            c.resize(len(ind), 1)

            # Weird stuff here because I am not familiar enough with numpy to make this cleaner
            # Basically need to take subset of numpy array but may not be sequential (ie rows 2 and 4 so can't use ':' in slicing)

            ind_pairs = list(itertools.product(ind.tolist(), repeat=2))
            nind_pairs = list(itertools.product(nind.tolist(), repeat=2))
            ind_pairs_list = []
            nind_pairs_list = []
            ind_nind_pairs_list = []
            for h in ind_pairs:
                ind_pairs_list.append(list(h))
            for h in nind_pairs:
                nind_pairs_list.append(list(h))
            for h in ind:
                for j in nind:
                    ind_nind_pairs_list.append([h, j])

            if len(ind) > 0:
                x, y = np.transpose(np.array(ind_pairs_list))
                D = C[x, y]
                D.resize((len(ind), len(ind)))
            else:
                D = np.array([[]])

            if len(nind) > 0:
                a, b = np.transpose(np.array(nind_pairs_list))
                A = C[a, b]
                A.resize((len(nind), len(nind)))
            else:
                A = np.array([[]])

            if len(ind) > 0 and len(nind) > 0:
                s, t = np.transpose(np.array(ind_nind_pairs_list))
                B = C[s, t]
                B.resize((len(ind), len(nind)))
            elif len(ind) > 0 and len(nind) == 0:
                B = np.zeros((len(ind), 1))
            elif len(nind) > 0 and len(ind) == 0:
                B = np.zeros((1, len(nind)))
            else:
                B = np.array([[]])

            if D.size > 0:
                Dinv = np.linalg.inv(D)
                # NOTE THROUGHOUT - order of transposes are sometimes different to original code because numpy doesn't always preserve 'rows' and 'columns'
                # the same as Octave does....
                implied_variance = np.matmul(np.matmul(np.transpose(c), Dinv), c)
                if implied_variance > d * p:
                    max_cov_value2 = np.diag(D) * d
                    c.resize(1, c.size)
                    replace_idx = (c**2 >= max_cov_value2 * p).nonzero()
                    if len(replace_idx) > 0:
                        max_cov_value2 = np.array([max_cov_value2])
                        c[replace_idx] = np.sign(c[replace_idx]) * np.sqrt(
                            max_cov_value2[replace_idx]
                        )
                        implied_variance = np.matmul(
                            np.matmul(c, Dinv), np.transpose(c)
                        )
                    if implied_variance > d * p:
                        c = np.sqrt(d * p / implied_variance) * c

                try:
                    x = np.matmul(np.matmul(np.transpose(B), Dinv), np.transpose(c))
                except ValueError:
                    x = np.matmul(np.matmul(np.transpose(B), Dinv), c)
            else:
                x = np.zeros((B.shape[1], 1))
            tmp = np.transpose(c)
            ctransp1D = tmp.reshape(-1)
            c1D = c.reshape(-1)
            tmp = np.transpose(x)
            xtransp1D = tmp.reshape(-1)
            x1D = x.reshape(-1)

            C[ind, i] = ctransp1D
            C[i, ind] = c1D
            if len(x) > 0:
                C[i, nind] = x1D
                C[nind, i] = xtransp1D

        return C

    def weighted_cchalf(
        self, this, other, assume_index_matching=False, use_binning=False, weighted=True
    ):
        if not use_binning:
            assert other.indices().size() == this.indices().size()
            if this.data().size() == 0:
                return None, None

            if assume_index_matching:
                (o, c) = (this, other)
            else:
                (o, c) = this.common_sets(other=other, assert_no_singles=True)

            # The case where the denominator is less or equal to zero is
            # pathological and should never arise in practice.
            if weighted:
                assert len(o.sigmas())
                assert len(c.sigmas())
                n = len(o.data())
                if n == 1:
                    return None, 1
                v_o = o.sigmas() ** 2
                v_c = c.sigmas() ** 2
                var_w = v_o + v_c
                joint_w = 1.0 / var_w
                sumjw = flex.sum(joint_w)
                norm_jw = joint_w / sumjw
                xbar = flex.sum(o.data() * norm_jw)
                ybar = flex.sum(c.data() * norm_jw)
                sxy = flex.sum((o.data() - xbar) * (c.data() - ybar) * norm_jw)

                sx = flex.sum((o.data() - xbar) ** 2 * norm_jw)
                sy = flex.sum((c.data() - ybar) ** 2 * norm_jw)
                # what is neff? neff = 1/V2
                # V2 = flex.sum(norm_jw**2)
                # use entropy based approach
                neff = math.exp(-1.0 * flex.sum(norm_jw * flex.log(norm_jw)))
                # print(n, neff, 1.0/V2)
                return (sxy / ((sx * sy) ** 0.5), neff)
            else:
                n = len(o.data())
                xbar = flex.sum(o.data()) / n
                ybar = flex.sum(c.data()) / n
                sxy = flex.sum((o.data() - xbar) * (c.data() - ybar))
                sx = flex.sum((o.data() - xbar) ** 2)
                sy = flex.sum((c.data() - ybar) ** 2)

                return (sxy / ((sx * sy) ** 0.5), n)
        assert this.binner() is not None
        results = []
        n_eff = []
        for i_bin in this.binner().range_all():
            sel = this.binner().selection(i_bin)
            cchalf, neff = self.weighted_cchalf(
                this.select(sel),
                other.select(sel),
                assume_index_matching=assume_index_matching,
                use_binning=False,
                weighted=weighted,
            )
            results.append(cchalf)
            n_eff.append(neff)
        return binned_data(binner=this.binner(), data=results, data_fmt="%7.4f"), n_eff
