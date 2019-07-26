"""Methods for symmetry determination from partial datasets.

This module implements the methods of `Gildea, R. J. & Winter, G. (2018).
Acta Cryst. D74, 405-410 <https://doi.org/10.1107/S2059798318002978>`_ for
determination of Patterson group symmetry from sparse multi-crystal data sets in
the presence of an indexing ambiguity.
"""
from __future__ import absolute_import, division, print_function

import copy
import json
import logging
import math
from collections import OrderedDict

import iotbx.phil
from cctbx import sgtbx
from dials.algorithms.indexing.symmetry import find_matching_symmetry
from dials.algorithms.symmetry.cosym import target
from dials.algorithms.symmetry.cosym import engine
from dials.algorithms.symmetry import symmetry_base
from dials.algorithms.symmetry.determine_space_group import ScoreCorrelationCoefficient
from dials.util.observer import Subject
from libtbx import Auto
from libtbx import table_utils
from scitbx import matrix
from scitbx.array_family import flex

logger = logging.getLogger(__name__)

phil_scope = iotbx.phil.parse(
    """\

normalisation = kernel quasi *ml_iso ml_aniso
  .type = choice

d_min = Auto
  .type = float(value_min=0)

min_i_mean_over_sigma_mean = 4
  .type = float(value_min=0)

min_cc_half = 0.6
  .type = float(value_min=0, value_max=1)

lattice_group = None
  .type = space_group

space_group = None
  .type = space_group

dimensions = Auto
  .type = int(value_min=2)

use_curvatures = True
  .type = bool

weights = count standard_error
  .type = choice

min_pairs = 3
  .type = int(value_min=1)
  .help = 'Minimum number of pairs for inclusion of correlation coefficient in calculation of Rij matrix.'

termination_params {
  max_iterations = 100
    .type = int(value_min=0)
  max_calls = None
    .type = int(value_min=0)
  traditional_convergence_test = True
    .type = bool
  traditional_convergence_test_eps = 1
    .type = float
  drop_convergence_test_n_test_points=5
    .type = int(value_min=2)
  drop_convergence_test_max_drop_eps=1.e-5
    .type = float(value_min=0)
  drop_convergence_test_iteration_coefficient=2
    .type = float(value_min=1)
}

cluster {
  method = dbscan bisect minimize_divide agglomerative *seed
    .type = choice
  n_clusters = auto
    .type = int(value_min=1)
  dbscan {
    eps = 0.5
      .type = float(value_min=0)
    min_samples = 5
      .type = int(value_min=1)
  }
  bisect {
    axis = 0
      .type = int(value_min=0)
  }
  seed {
    min_silhouette_score = 0.2
      .type = float(value_min=-1, value_max=1)
  }
}

nproc = 1
  .type = int(value_min=1)
  .help = "The number of processes to use."

"""
)


class CosymAnalysis(symmetry_base, Subject):
    """Peform cosym analysis.

    Peform cosym analysis on the input intensities using the methods of
    `Gildea, R. J. & Winter, G. (2018). Acta Cryst. D74, 405-410
    <https://doi.org/10.1107/S2059798318002978>`_ for
    determination of Patterson group symmetry from sparse multi-crystal data sets in
    the presence of an indexing ambiguity.

    """

    def __init__(self, intensities, params):
        """Initialise a CosymAnalysis object.

        Args:
          intensities (cctbx.miller.array): The intensities on which to perform
            cosym anaylsis.
          params (libtbx.phil.scope_extract): Parameters for the analysis.

        """
        super(CosymAnalysis, self).__init__(
            intensities,
            normalisation=params.normalisation,
            lattice_symmetry_max_delta=5.0,
            d_min=params.d_min,
            min_i_mean_over_sigma_mean=params.min_i_mean_over_sigma_mean,
            min_cc_half=params.min_cc_half,
            relative_length_tolerance=None,
            absolute_angle_tolerance=None,
        )
        Subject.__init__(
            self, events=["optimised", "analysed_symmetry", "analysed_clusters"]
        )

        self.params = params
        if self.params.space_group is not None:

            def _map_space_group_to_input_cell(intensities, space_group):
                best_subgroup = find_matching_symmetry(
                    intensities.unit_cell(), space_group
                )
                cb_op_inp_best = best_subgroup["cb_op_inp_best"]
                best_subsym = best_subgroup["best_subsym"]
                cb_op_best_ref = best_subsym.change_of_basis_op_to_reference_setting()
                ref_subsym = best_subsym.change_basis(cb_op_best_ref)
                cb_op_ref_primitive = (
                    ref_subsym.change_of_basis_op_to_primitive_setting()
                )
                sg_cb_op_inp_primitive = (
                    space_group.info().change_of_basis_op_to_primitive_setting()
                )
                sg_primitive = space_group.change_basis(sg_cb_op_inp_primitive)
                sg_best = sg_primitive.change_basis(
                    (cb_op_ref_primitive * cb_op_best_ref).inverse()
                )
                # best_subgroup above is the bravais type, so create thin copy here with the
                # user-input space group instead
                best_subgroup = {
                    "best_subsym": best_subsym.customized_copy(
                        space_group_info=sg_best.info()
                    ),
                    "cb_op_inp_best": cb_op_inp_best,
                }

                intensities = intensities.customized_copy(
                    space_group_info=sg_best.change_basis(
                        cb_op_inp_best.inverse()
                    ).info()
                )
                return intensities, best_subgroup

            self.intensities, self.best_subgroup = _map_space_group_to_input_cell(
                self.intensities, self.params.space_group.group()
            )
            self.input_space_group = self.intensities.space_group()

        else:
            self.input_space_group = None

        if self.params.lattice_group is not None:
            tmp_intensities, _ = _map_space_group_to_input_cell(
                self.intensities, self.params.lattice_group.group()
            )
            self.params.lattice_group = tmp_intensities.space_group_info()

    def _intialise_target(self):
        if self.params.dimensions is Auto:
            dimensions = None
        else:
            dimensions = self.params.dimensions
        if self.params.lattice_group is not None:
            self.lattice_group = (
                self.params.lattice_group.group()
                .build_derived_patterson_group()
                .info()
                .primitive_setting()
                .group()
            )
        self.target = target.Target(
            self.intensities,
            self.dataset_ids,
            min_pairs=self.params.min_pairs,
            lattice_group=self.lattice_group,
            dimensions=dimensions,
            weights=self.params.weights,
            nproc=self.params.nproc,
        )

    def _determine_dimensions(self):
        if self.params.dimensions is Auto and self.target.dim == 2:
            self.params.dimensions = 2
        elif self.params.dimensions is Auto:
            dimensions = []
            functional = []
            explained_variance = []
            explained_variance_ratio = []
            for dim in range(1, self.target.dim + 1):
                self.target.set_dimensions(dim)
                self._optimise()
                logger.info("Functional: %g" % self.minimizer.f)
                self._principal_component_analysis()
                dimensions.append(dim)
                functional.append(self.minimizer.f)
                explained_variance.append(self.explained_variance)
                explained_variance_ratio.append(self.explained_variance_ratio)

            # Find the elbow point of the curve, in the same manner as that used by
            # distl spotfinder for resolution method 1 (Zhang et al 2006).
            # See also dials/algorithms/spot_finding/per_image_analysis.py

            x = flex.double(dimensions)
            y = flex.double(functional)
            slopes = (y[-1] - y[:-1]) / (x[-1] - x[:-1])
            p_m = flex.min_index(slopes)

            x1 = matrix.col((x[p_m], y[p_m]))
            x2 = matrix.col((x[-1], y[-1]))

            gaps = flex.double()
            v = matrix.col(((x2[1] - x1[1]), -(x2[0] - x1[0]))).normalize()

            for i in range(p_m, len(x)):
                x0 = matrix.col((x[i], y[i]))
                r = x1 - x0
                g = abs(v.dot(r))
                gaps.append(g)

            p_g = flex.max_index(gaps)

            x_g = x[p_g + p_m]

            logger.info("Best number of dimensions: %i" % x_g)
            self.target.set_dimensions(int(x_g))

    def run(self):
        self._intialise_target()
        self._determine_dimensions()
        self._optimise()
        self._principal_component_analysis()

        self._analyse_symmetry()
        self._cluster_analysis()

    @Subject.notify_event(event="optimised")
    def _optimise(self):
        NN = len(self.input_intensities)
        dim = self.target.dim
        n_sym_ops = len(self.target.get_sym_ops())
        coords = flex.random_double(NN * n_sym_ops * dim)

        import scitbx.lbfgs

        tp = self.params.termination_params
        termination_params = scitbx.lbfgs.termination_parameters(
            traditional_convergence_test=tp.traditional_convergence_test,
            traditional_convergence_test_eps=tp.traditional_convergence_test_eps,
            drop_convergence_test_n_test_points=tp.drop_convergence_test_n_test_points,
            drop_convergence_test_max_drop_eps=tp.drop_convergence_test_max_drop_eps,
            drop_convergence_test_iteration_coefficient=tp.drop_convergence_test_iteration_coefficient,
            # min_iterations=tp.min_iterations,
            max_iterations=tp.max_iterations,
            max_calls=tp.max_calls,
        )

        M = engine.lbfgs_with_curvs(
            self.target,
            coords,
            use_curvatures=self.params.use_curvatures,
            termination_params=termination_params,
        )
        self.minimizer = M

        coords = M.x.deep_copy()
        coords.reshape(flex.grid(dim, NN * n_sym_ops))
        coords.matrix_transpose_in_place()
        self.coords = coords

    def _principal_component_analysis(self):
        # Perform PCA
        from sklearn.decomposition import PCA

        X = self.coords.as_numpy_array()
        pca = PCA().fit(X)
        logger.info("Principal component analysis:")
        logger.info(
            "Explained variance: "
            + ", ".join(["%.2g" % v for v in pca.explained_variance_])
        )
        logger.info(
            "Explained variance ratio: "
            + ", ".join(["%.2g" % v for v in pca.explained_variance_ratio_])
        )
        self.explained_variance = pca.explained_variance_
        self.explained_variance_ratio = pca.explained_variance_ratio_
        if self.target.dim > 3:
            pca.n_components = 3
        x_reduced = pca.fit_transform(X)

        import numpy

        self.coords_reduced = flex.double(numpy.ascontiguousarray(x_reduced))

    @Subject.notify_event(event="analysed_symmetry")
    def _analyse_symmetry(self):
        if self.input_space_group is not None:
            self.best_solution = None
            self._symmetry_analysis = None
            return

        sym_ops = [
            sgtbx.rt_mx(s).new_denominators(1, 12) for s in self.target.get_sym_ops()
        ]
        self._symmetry_analysis = SymmetryAnalysis(
            self.coords, sym_ops, self.subgroups, self.cb_op_inp_min
        )
        logger.info(str(self._symmetry_analysis))
        self.best_solution = self._symmetry_analysis.best_solution
        self.best_subgroup = self.best_solution.subgroup

        cosets = sgtbx.cosets.left_decomposition(
            self.lattice_group, self.best_solution.subgroup["subsym"].space_group()
        )
        self.params.cluster.n_clusters = len(cosets.partitions)

    def _space_group_for_dataset(self, dataset_id, sym_ops):
        if self.input_space_group is not None:
            sg = copy.deepcopy(self.input_space_group)
        else:
            sg = sgtbx.space_group()
        ref_sym_op_id = None
        ref_cluster_id = None
        for sym_op_id in range(len(sym_ops)):
            i_cluster = self.cluster_labels[
                len(self.input_intensities) * sym_op_id + dataset_id
            ]
            if i_cluster < 0:
                continue
            if ref_sym_op_id is None:
                ref_sym_op_id = sym_op_id
                ref_cluster_id = i_cluster
                continue
            op = sym_ops[ref_sym_op_id].inverse().multiply(sym_ops[sym_op_id])
            if i_cluster == ref_cluster_id:
                sg.expand_smx(op.new_denominators(1, 12))
        return sg.make_tidy()

    def _reindexing_ops_for_dataset(self, dataset_id, sym_ops, cosets):
        reindexing_ops = {}
        # Number of clusters in labels, ignoring noise if present.
        n_clusters = len(set(self.cluster_labels)) - (
            1 if -1 in self.cluster_labels else 0
        )

        for i_cluster in range(n_clusters):
            isel = (self.cluster_labels == i_cluster).iselection()
            dataset_ids = isel % len(self.input_intensities)
            sel = (dataset_ids == dataset_id).iselection()
            for s in sel:
                sym_op_id = isel[s] // len(self.input_intensities)
                for partition in cosets.partitions:
                    if sym_ops[sym_op_id] in partition:
                        if i_cluster not in reindexing_ops:
                            cb_op = sgtbx.change_of_basis_op(
                                partition[0]
                            ).new_denominators(self.cb_op_inp_min)
                            reindexing_ops[i_cluster] = (
                                self.cb_op_inp_min.inverse()
                                * cb_op
                                * self.cb_op_inp_min
                            ).as_xyz()

        return reindexing_ops

    @Subject.notify_event(event="analysed_clusters")
    def _cluster_analysis(self):

        if self.params.cluster.n_clusters == 1:
            self.cluster_labels = flex.double(self.coords.all()[0])
        else:
            self.cluster_labels = self._do_clustering(self.params.cluster.method)

        sym_ops = [
            sgtbx.rt_mx(s).new_denominators(1, 12) for s in self.target.get_sym_ops()
        ]

        reindexing_ops = {}
        space_groups = {}

        for dataset_id in range(len(self.input_intensities)):
            space_groups[dataset_id] = self._space_group_for_dataset(
                dataset_id, sym_ops
            )

            cosets = sgtbx.cosets.left_decomposition(
                self.target._lattice_group, space_groups[dataset_id]
            )

            reindexing_ops[dataset_id] = self._reindexing_ops_for_dataset(
                dataset_id, sym_ops, cosets
            )

        self.space_groups = space_groups
        self.reindexing_ops = reindexing_ops

    def _do_clustering(self, method):
        if method == "dbscan":
            clustering = self._dbscan_clustering
        elif method == "bisect":
            clustering = self._bisect_clustering
        elif method == "minimize_divide":
            clustering = self._minimize_divide_clustering
        elif method == "agglomerative":
            clustering = self._agglomerative_clustering
        elif method == "seed":
            clustering = self._seed_clustering
        return clustering()

    def _dbscan_clustering(self):
        from sklearn.preprocessing import StandardScaler

        X = self.coords_reduced.as_numpy_array()
        X = StandardScaler().fit_transform(X)

        # Perform cluster analysis
        from sklearn.cluster import DBSCAN

        db = DBSCAN(
            eps=self.params.cluster.dbscan.eps,
            min_samples=self.params.cluster.dbscan.min_samples,
        ).fit(X)
        import numpy as np

        return flex.int(db.labels_.astype(np.int32))

    def _bisect_clustering(self):
        assert self.params.cluster.n_clusters in (2, Auto)
        axis = self.params.cluster.bisect.axis
        assert axis < self.coords_reduced.all()[1]
        x = self.coords_reduced[:, axis : axis + 1].as_1d()
        cluster_labels = flex.int(x.size(), 0)
        cluster_labels.set_selected(x > 0, 1)
        return cluster_labels

    def _minimize_divide_clustering(self):
        assert self.params.cluster.n_clusters in (2, Auto)
        x = self.coords_reduced[:, :1].as_1d()
        y = self.coords_reduced[:, 1:2].as_1d()
        from cctbx.merging.brehm_diederichs import minimize_divide

        selection = minimize_divide(x, y).plus_minus()
        cluster_labels = flex.int(x.size(), 0)
        cluster_labels.set_selected(selection, 1)
        return cluster_labels

    def _agglomerative_clustering(self):
        X = self.coords.as_numpy_array()

        # Perform cluster analysis
        from sklearn.cluster import AgglomerativeClustering
        import numpy as np

        model = AgglomerativeClustering(
            n_clusters=self.params.cluster.n_clusters,
            linkage="average",
            affinity="cosine",
        )
        model.fit(X)
        return flex.int(model.labels_.astype(np.int32))

    def _seed_clustering(self):
        from dials.algorithms.symmetry.cosym.seed_clustering import seed_clustering

        clustering = seed_clustering(
            self.coords,
            len(self.input_intensities),
            len(self.target.get_sym_ops()),
            min_silhouette_score=self.params.cluster.seed.min_silhouette_score,
            n_clusters=self.params.cluster.n_clusters,
        )
        return clustering.cluster_labels

    def as_dict(self):
        """Return a dictionary representation of the results.

        Returns:
          dict

        """
        d = {
            "input_symmetry": {
                "hall_symbol": self.input_intensities[0]
                .space_group()
                .type()
                .hall_symbol(),
                "unit_cell": self.median_unit_cell.parameters(),
            },
            "cb_op_inp_min": self.cb_op_inp_min.as_xyz(),
            "min_cell_symmetry": {
                "hall_symbol": self.intensities.space_group().type().hall_symbol(),
                "unit_cell": self.intensities.unit_cell().parameters(),
            },
            "lattice_point_group": self.lattice_group.type().hall_symbol(),
        }

        if self._symmetry_analysis is not None:
            d.update(self._symmetry_analysis.as_dict())
        return d

    def as_json(self, filename=None, indent=2):
        """Return a json representation of the results.

        Args:
          filename (str): Optional filename to export the json representation of
            the results.
          indent (int): The indent level for pretty-printing of the json. If ``None``
            is the most compact representation.

        Returns:
          str:

        """
        d = self.as_dict()

        json_str = json.dumps(d, indent=indent)
        if filename:
            with open(filename, "w") as f:
                f.write(json_str)
        return json.dumps(d, indent=indent)


class SymmetryAnalysis(object):
    def __init__(self, coords, sym_ops, subgroups, cb_op_inp_min):

        import scipy.spatial.distance as ssd

        self.subgroups = subgroups
        self.cb_op_inp_min = cb_op_inp_min
        X = coords.as_numpy_array()
        n_datasets = coords.all()[0] // len(sym_ops)
        dist_mat = ssd.pdist(X, metric="cosine")
        cos_angle = 1 - ssd.squareform(dist_mat)

        self._sym_ops_cos_angle = OrderedDict()
        for dataset_id in range(n_datasets):
            for ref_sym_op_id in range(len(sym_ops)):
                ref_idx = n_datasets * ref_sym_op_id + dataset_id
                for sym_op_id in range(ref_sym_op_id + 1, len(sym_ops)):
                    op = sym_ops[ref_sym_op_id].inverse().multiply(sym_ops[sym_op_id])
                    op = op.new_denominators(1, 12)
                    comp_idx = n_datasets * sym_op_id + dataset_id
                    self._sym_ops_cos_angle.setdefault(op, flex.double())
                    self._sym_ops_cos_angle[op].append(cos_angle[ref_idx, comp_idx])

        self._score_symmetry_elements()
        self._score_laue_groups()

    def _score_symmetry_elements(self):
        self.sym_op_scores = OrderedDict()
        for op, cos_angle in self._sym_ops_cos_angle.items():
            cc_true = 1
            cc = flex.mean(cos_angle)
            score = ScoreSymmetryElement(cc, sigma_cc=0.1, cc_true=cc_true)
            score.sym_op = op
            self.sym_op_scores[op] = score

    def _score_laue_groups(self):
        subgroup_scores = [
            ScoreSubGroup(subgrp, list(self.sym_op_scores.values()))
            for subgrp in self.subgroups.result_groups
        ]
        total_likelihood = sum(score.likelihood for score in subgroup_scores)
        for score in subgroup_scores:
            score.likelihood /= total_likelihood
        self.subgroup_scores = sorted(
            subgroup_scores, key=lambda score: score.likelihood, reverse=True
        )

        # The 'confidence' scores are derived from the total probability of the best
        # solution p_best and that for the next best solution p_next:
        #   confidence = [p_best * (p_best - p_next)]^1/2.

        for i, score in enumerate(self.subgroup_scores[:-1]):
            next_score = self.subgroup_scores[i + 1]
            if score.likelihood > 0 and next_score.likelihood > 0:
                lgc = score.likelihood * (score.likelihood - next_score.likelihood)
                confidence = abs(lgc) ** 0.5
                if lgc < 0:
                    confidence = -confidence
                score.confidence = confidence

        self.best_solution = self.subgroup_scores[0]

    @staticmethod
    def sym_ops_table(d):
        header = ("likelihood", "Z-CC", "CC", "", "Operator")
        rows = [header]
        for score in d["sym_op_scores"]:
            rows.append(
                (
                    "%.3f" % score["likelihood"],
                    "%.2f" % score["z_cc"],
                    "%.2f" % score["cc"],
                    score["stars"],
                    str(sgtbx.rt_mx(str(score["operator"])).r().info()),
                )
            )
        return rows

    @staticmethod
    def subgroups_table(d):
        header = (
            "Patterson group",
            "",
            "Likelihood",
            "NetZcc",
            "Zcc+",
            "Zcc-",
            "delta",
            "Reindex operator",
        )
        rows = [header]
        for score in d["subgroup_scores"]:
            rows.append(
                (
                    str(
                        sgtbx.space_group(
                            hall_symbol=str(score["patterson_group"])
                        ).info()
                    ),
                    score["stars"],
                    "%.3f" % score["likelihood"],
                    "% .2f" % score["z_cc_net"],
                    "% .2f" % score["z_cc_for"],
                    "% .2f" % score["z_cc_against"],
                    "%.1f" % score["max_angular_difference"],
                    str(sgtbx.change_of_basis_op(str(score["cb_op"]))),
                )
            )
        return rows

    @staticmethod
    def summary_table(d):
        best_subgroup = d["subgroup_scores"][0]
        return (
            (
                "Best solution",
                str(
                    sgtbx.space_group(
                        hall_symbol=str(best_subgroup["patterson_group"])
                    ).info()
                ),
            ),
            (
                "Unit cell",
                "%.3f %.3f %.3f %.1f %.1f %.1f" % tuple(best_subgroup["unit_cell"]),
            ),
            ("Reindex operator", best_subgroup["cb_op"]),
            ("Laue group probability", "%.3f" % best_subgroup["likelihood"]),
            ("Laue group confidence", "%.3f" % best_subgroup["confidence"]),
        )

    def __str__(self):
        """Return a string representation of the results.

        Returns:
          str:

        """
        output = []
        output.append("Scoring individual symmetry elements")
        d = self.as_dict()
        output.append(
            table_utils.format(self.sym_ops_table(d), has_header=True, delim="  ")
        )

        output.append("Scoring all possible sub-groups")
        output.append(
            table_utils.format(self.subgroups_table(d), has_header=True, delim="  ")
        )

        output.append(
            "Best solution: %s"
            % self.best_solution.subgroup["best_subsym"].space_group_info()
        )
        output.append(
            "Unit cell: %s" % self.best_solution.subgroup["best_subsym"].unit_cell()
        )
        output.append(
            "Reindex operator: %s"
            % (self.best_solution.subgroup["cb_op_inp_best"] * self.cb_op_inp_min)
        )
        output.append("Laue group probability: %.3f" % self.best_solution.likelihood)
        output.append("Laue group confidence: %.3f" % self.best_solution.confidence)
        return "\n".join(output)

    def as_dict(self):
        """Return a dictionary representation of the results.

        Returns:
          dict

        """
        d = {"cb_op_inp_min": self.cb_op_inp_min.as_xyz()}

        d["sym_op_scores"] = []
        for rt_mx, score in self.sym_op_scores.items():
            dd = score.as_dict()
            dd["operator"] = rt_mx.as_xyz()
            d["sym_op_scores"].append(dd)

        d["subgroup_scores"] = []
        for score in self.subgroup_scores:
            dd = score.as_dict()
            dd["cb_op"] = (
                sgtbx.change_of_basis_op(dd["cb_op"]) * self.cb_op_inp_min
            ).as_xyz()
            d["subgroup_scores"].append(dd)

        return d


class ScoreSymmetryElement(object):
    """Analyse intensities for presence of a given symmetry operation.

    1) Calculate the probability of observing this CC if the sym op is present,
       p(CC; S), modelled by a Cauchy distribution centred on cc_true and width
       gamma = sigma_cc.

    2) Calculate the probability of observing this CC if the sym op is
       NOT present, p(CC; !S).

    3) Calculate the likelihood of symmetry element being present,
       p(S; CC) = p(CC; S) / (p(CC; S) + p(CC; !S))

    See appendix A1 of `Evans, P. R. (2011). Acta Cryst. D67, 282-292.
    <https://doi.org/10.1107/S090744491003982X>`_

    """

    def __init__(self, cc, sigma_cc, cc_true):
        """Initialise a ScoreSymmetryElement object.

        Args:
          cc (float): the correlation coefficient for this symmetry element
          sigma_cc (float): the estimated error in the correlation coefficient
          cc_true (float): the expected value of CC if the symmetry element is present,
            E(CC; S)

        """

        self.cc = cc
        self.sigma_cc = sigma_cc
        self.z_cc = self.cc / self.sigma_cc
        score_cc = ScoreCorrelationCoefficient(self.cc, self.sigma_cc, cc_true)
        self.p_cc_given_s = score_cc.p_cc_given_s
        self.p_cc_given_not_s = score_cc.p_cc_given_not_s
        self.likelihood = score_cc.p_s_given_cc

    @property
    def stars(self):
        # define stars attribute - used mainly for output
        if self.likelihood > 0.9:
            stars = "***"
        elif self.likelihood > 0.7:
            stars = "**"
        elif self.likelihood > 0.5:
            stars = "*"
        else:
            stars = ""
        return stars

    def as_dict(self):
        """Return a dictionary representation of the symmetry element scoring.

        The dictionary will contain the following keys:
          - likelihood: The likelihood of the symmetry element being present
          - z_cc: The Z-score for the correlation coefficent
          - cc: The correlation coefficient for the symmetry element
          - operator: The xyz representation of the symmetry element

        Returns:
          dict:

        """

        return {
            "likelihood": self.likelihood,
            "z_cc": self.z_cc,
            "cc": self.cc,
            "stars": self.stars,
        }


class ScoreSubGroup(object):
    """Score the probability of a given subgroup being the true subgroup.

    1) Calculates overall Zcc scores for symmetry elements present/absent from
       the subgroup.

    2) Calculates the overall likelihood for this subgroup.

    See appendix A2 of `Evans, P. R. (2011). Acta Cryst. D67, 282-292.
    <https://doi.org/10.1107/S090744491003982X>`_

    """

    def __init__(self, subgroup, sym_op_scores):
        """Initialise a ScoreSubGroup object.

        Args:
          subgroup (dict): A dictionary describing the subgroup as generated by
            :class:`cctbx.sgtbx.lattice_symmetry.metric_subgroups`.
          sym_op_scores (list): A list of :class:`ScoreSymmetryElement` objects for each
            symmetry element possibly in the lattice symmetry.

        """
        # Combined correlation coefficients for symmetry operations
        # present/absent from subgroup
        self.subgroup = subgroup
        patterson_group = subgroup["subsym"].space_group()

        # Overall Zcc scores for symmetry elements present/absent from subgroup
        self.z_cc_for = 0
        self.z_cc_against = 0
        n_for = 0
        n_against = 0
        PL_for = 0
        PL_against = 0
        power = 2
        for score in sym_op_scores:
            if score.sym_op in patterson_group:
                self.z_cc_for += score.z_cc ** power
                n_for += 1
                PL_for += math.log(score.p_cc_given_s)
            else:
                self.z_cc_against += score.z_cc ** power
                n_against += 1
                PL_against += math.log(score.p_cc_given_not_s)

        # Overall likelihood for this subgroup
        self.likelihood = math.exp(PL_for + PL_against)

        if n_against > 0:
            self.z_cc_against = (self.z_cc_against / n_against) ** (1 / power)
        if n_for > 0:
            self.z_cc_for = (self.z_cc_for / n_for) ** (1 / power)
        self.z_cc_net = self.z_cc_for - self.z_cc_against
        self.confidence = 0

    def __str__(self):
        """Return a string representation of the subgroup scores.

        Returns:
          str:

        """
        return "%s %.3f %.2f %.2f %.2f" % (
            self.subgroup["best_subsym"].space_group_info(),
            self.likelihood,
            self.z_cc_net,
            self.z_cc_for,
            self.z_cc_against,
        )

    @property
    def stars(self):
        if self.likelihood > 0.8:
            stars = "***"
        elif self.likelihood > 0.6:
            stars = "**"
        elif self.likelihood > 0.4:
            stars = "*"
        else:
            stars = ""
        return stars

    def as_dict(self):
        """Return a dictionary representation of the subgroup scoring.

        The dictionary will contain the following keys:
          - patterson_group: The current subgroup
          - likelihood: The likelihood of the subgroup being correct
          - confidence: The confidence of the subgroup being correct
          - z_cc_for: The combined Z-scores for all symmetry elements present in the
            subgroup
          - z_cc_against: The combined Z-scores for all symmetry elements present in
            the lattice group but not in the subgroup
          - z_cc_net: The net Z-score, i.e. z_cc_for - z_cc_against
          - max_angular_difference: The maximum angular difference between the
            symmetrised unit cell and the P1 unit cell.
          - cb_op: The change of basis operation from the input unit cell to the
            'best' unit cell.

        Returns:
          dict:

        """
        return {
            "patterson_group": self.subgroup["best_subsym"]
            .space_group()
            .type()
            .hall_symbol(),
            "unit_cell": self.subgroup["best_subsym"].unit_cell().parameters(),
            "likelihood": self.likelihood,
            "confidence": self.confidence,
            "z_cc_net": self.z_cc_net,
            "z_cc_for": self.z_cc_for,
            "z_cc_against": self.z_cc_against,
            "max_angular_difference": self.subgroup["max_angular_difference"],
            "cb_op": "%s" % (self.subgroup["cb_op_inp_best"]),
            "stars": self.stars,
        }
