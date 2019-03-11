"""Methods for symmetry determination from partial datasets.

This module implements the methods of `Gildea, R. J. & Winter, G. (2018).
Acta Cryst. D74, 405-410 <https://doi.org/10.1107/S2059798318002978>`_ for
determination of Patterson group symmetry from sparse multi-crystal data sets in
the presence of an indexing ambiguity.
"""
from __future__ import absolute_import, division, print_function

import logging

logger = logging.getLogger(__name__)

import copy
from collections import OrderedDict
import math

from libtbx import Auto
from libtbx import table_utils
from scitbx.array_family import flex
from scitbx import matrix
from cctbx import sgtbx
import iotbx.phil

from dials.algorithms.symmetry.cosym import target
from dials.algorithms.symmetry.cosym import engine
from dials.algorithms.symmetry import symmetry_base

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

dimensions = None
  .type = int(value_min=2)

use_curvatures = True
  .type = bool

weights = count standard_error
  .type = choice

min_pairs = 3
  .type = int(value_min=1)
  .help = 'Minimum number of pairs for inclusion of correlation coefficient in calculation of Rij matrix.'

save_plot = True
  .type = bool

plot_prefix = ''
  .type = str

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


class analyse_datasets(symmetry_base):
    """Peform cosym analysis.

    Peform cosym analysis on the input intensities using the methods of
    `Gildea, R. J. & Winter, G. (2018). Acta Cryst. D74, 405-410
    <https://doi.org/10.1107/S2059798318002978>`_ for
    determination of Patterson group symmetry from sparse multi-crystal data sets in
    the presence of an indexing ambiguity.

    """

    def __init__(self, intensities, params):
        """Initialise an analyse_datasets object.

        Args:
          intensities (cctbx.miller.array): The intensities on which to perform
            cosym anaylsis.
          params (libtbx.phil.scope_extract): Parameters for the analysis.

        """
        self.input_space_group = intensities[0].space_group()
        super(analyse_datasets, self).__init__(
            intensities,
            normalisation=params.normalisation,
            lattice_symmetry_max_delta=5.0,
            d_min=params.d_min,
            min_i_mean_over_sigma_mean=params.min_i_mean_over_sigma_mean,
            min_cc_half=params.min_cc_half,
            relative_length_tolerance=None,
            absolute_angle_tolerance=None,
        )

        self.params = params
        self.intensities = self.intensities.customized_copy(
            space_group_info=self.input_space_group.change_basis(
                self.cb_op_inp_min
            ).info()
        )
        if self.params.dimensions is Auto:
            dimensions = None
        else:
            dimensions = self.params.dimensions
        lattice_group = None
        if self.params.lattice_group is not None:
            lattice_group = (
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
            lattice_group=lattice_group,
            dimensions=dimensions,
            weights=self.params.weights,
            nproc=self.params.nproc,
        )
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

            p_k = flex.max_index(gaps)
            g_k = gaps[p_k]
            p_g = p_k

            x_g = x[p_g + p_m]
            y_g = y[p_g + p_m]

            logger.info("Best number of dimensions: %i" % x_g)
            self.target.set_dimensions(int(x_g))

            if params.save_plot:
                from matplotlib import pyplot as plt

                fig = plt.figure(figsize=(10, 8))
                plt.clf()
                plt.plot(dimensions, functional)
                plt.plot([x_g, x_g], plt.ylim())
                plt.xlabel("Dimensions")
                plt.ylabel("Functional")
                plt.savefig("%sfunctional_vs_dimension.png" % params.plot_prefix)

                plt.clf()
                for dim, expl_var in zip(dimensions, explained_variance):
                    plt.plot(range(1, dim + 1), expl_var, label="%s" % dim)
                plt.plot([x_g, x_g], plt.ylim())
                plt.xlabel("Dimension")
                plt.ylabel("Explained variance")
                plt.savefig(
                    "%sexplained_variance_vs_dimension.png" % params.plot_prefix
                )

                plt.clf()
                for dim, expl_var_ratio in zip(dimensions, explained_variance_ratio):
                    plt.plot(range(1, dim + 1), expl_var_ratio, label="%s" % dim)
                plt.plot([x_g, x_g], plt.ylim())
                plt.xlabel("Dimension")
                plt.ylabel("Explained variance ratio")
                plt.savefig(
                    "%sexplained_variance_ratio_vs_dimension.png" % params.plot_prefix
                )
                plt.close(fig)

        self._optimise()
        self._principal_component_analysis()

        self._analyse_symmetry()
        self._cluster_analysis()
        if self.params.save_plot:
            self._plot()

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

    def _analyse_symmetry(self):
        if self.input_space_group.type().number() > 1:
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

        cosets = sgtbx.cosets.left_decomposition(
            self.lattice_group, self.best_solution.subgroup["subsym"].space_group()
        )
        self.params.cluster.n_clusters = len(cosets.partitions)

    def _cosine_analysis(self):
        from scipy.cluster import hierarchy
        import scipy.spatial.distance as ssd

        X = self.coords.as_numpy_array()
        dist_mat = ssd.pdist(X, metric="cosine")
        cos_angle = 1 - ssd.squareform(dist_mat)
        linkage_matrix = hierarchy.linkage(dist_mat, method="average")

        c, coph_dists = hierarchy.cophenet(linkage_matrix, dist_mat)
        logger.debug(
            "Cophenetic correlation coefficient between heirarchical clustering and pairwise distance matrix: %.3f"
            % c
        )

        if self.params.save_plot:
            plot_matrix(
                cos_angle,
                linkage_matrix,
                "%scos_angle_matrix.png" % self.params.plot_prefix,
            )
            plot_dendrogram(
                linkage_matrix, "%scos_angle_dendrogram.png" % self.params.plot_prefix
            )

        sym_ops = [
            sgtbx.rt_mx(s).new_denominators(1, 12) for s in self.target.get_sym_ops()
        ]

        sym_ops_cos_angle = OrderedDict()

        for dataset_id in range(len(self.input_intensities)):
            ref_sym_op_id = None
            ref_cluster_id = None
            for sym_op_id in range(len(sym_ops)):
                if ref_sym_op_id is None:
                    ref_sym_op_id = sym_op_id
                    continue
                op = sym_ops[ref_sym_op_id].inverse().multiply(sym_ops[sym_op_id])
                op = op.new_denominators(1, 12)

                ref_idx = len(self.input_intensities) * ref_sym_op_id + dataset_id
                comp_idx = len(self.input_intensities) * sym_op_id + dataset_id
                sym_ops_cos_angle.setdefault(op, flex.double())
                sym_ops_cos_angle[op].append(cos_angle[ref_idx, comp_idx])

        # print symops sorted by average cos(angle)
        sg = copy.deepcopy(self.input_space_group)
        rows = [["symop", "order", "sg", "mean(cos(angle))", "median(cos(angle))"]]
        perm = flex.sort_permutation(
            flex.double([flex.mean(ca) for ca in sym_ops_cos_angle.values()]),
            reverse=True,
        )
        for p in perm:
            op, ca = sym_ops_cos_angle.items()[p]
            sg.expand_smx(op)
            rows.append(
                (
                    str(op),
                    str(op.r().order()),
                    str(sg.info().reference_setting()),
                    "%.3f" % flex.mean(ca),
                    "%.3f" % flex.median(ca),
                )
            )
        logger.info(
            "Analysis of cos(angle) between points corresponding to the same datasets:"
        )
        logger.info(table_utils.format(rows, has_header=True))

    def _space_group_for_dataset(self, dataset_id, sym_ops):
        sg = copy.deepcopy(self.input_space_group)
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
            idx = flex.first_index(dataset_ids, dataset_id)
            sel = (dataset_ids == dataset_id).iselection()
            if idx >= 0:
                sym_op_id = isel[idx] // len(self.input_intensities)
            for s in sel:
                sym_op_id = isel[s] // len(self.input_intensities)
                for partition in cosets.partitions:
                    if sym_ops[sym_op_id] in partition:
                        if i_cluster not in reindexing_ops:
                            cb_op = sgtbx.change_of_basis_op(
                                partition[0]
                            ).new_denominators(self.cb_op_inp_min)
                            reindexing_ops[i_cluster] = cb_op.as_xyz()

        return reindexing_ops

    def _cluster_analysis(self):

        if self.params.cluster.n_clusters == 1:
            self.cluster_labels = flex.double(self.coords.all()[0])
        else:
            self.cluster_labels = self._do_clustering(self.params.cluster.method)

        cluster_miller_arrays = []

        space_groups = []

        sym_ops = [
            sgtbx.rt_mx(s).new_denominators(1, 12) for s in self.target.get_sym_ops()
        ]
        self.space_groups = space_groups

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
            plot_prefix=self.params.plot_prefix if self.params.save_plot else None,
        )
        return clustering.cluster_labels

    def _plot(self):
        self.target.plot_rij_matrix(plot_name="%srij.png" % self.params.plot_prefix)
        self.target.plot_rij_histogram(
            plot_name="%srij_hist.png" % self.params.plot_prefix
        )
        self.target.plot_rij_cumulative_frequency(
            plot_name="%srij_sorted.png" % self.params.plot_prefix
        )
        self.target.plot_wij_matrix(plot_name="%swij.png" % self.params.plot_prefix)
        self.target.plot_wij_histogram(
            plot_name="%swij_hist.png" % self.params.plot_prefix
        )
        self.target.plot_wij_cumulative_frequency(
            plot_name="%swij_sorted.png" % self.params.plot_prefix
        )

        coord_x = self.coords[:, 0:1].as_1d()
        coord_y = self.coords[:, 1:2].as_1d()
        assert coord_x.size() == coord_y.size(), (coord_x.size(), coord_y.size())
        coord_reduced_x = self.coords_reduced[:, 0:1].as_1d()
        coord_reduced_y = self.coords_reduced[:, 1:2].as_1d()
        _plot(
            (coord_x, coord_y),
            labels=self.cluster_labels,
            plot_name="%sxy.png" % self.params.plot_prefix,
        )
        _plot(
            (coord_reduced_x, coord_reduced_y),
            labels=self.cluster_labels,
            plot_name="%sxy_pca.png" % self.params.plot_prefix,
        )
        _plot_angles(
            (coord_x, coord_y),
            labels=self.cluster_labels,
            plot_name="%sphi_r.png" % self.params.plot_prefix,
        )
        _plot_angles(
            (coord_reduced_x, coord_reduced_y),
            labels=self.cluster_labels,
            plot_name="%sphi_r_pca.png" % self.params.plot_prefix,
        )

        if self.coords_reduced.all()[1] > 2:
            coord_z = self.coords[:, 2:3].as_1d()
            coord_reduced_z = self.coords_reduced[:, 2:3].as_1d()
            _plot(
                (coord_x, coord_y, coord_z),
                labels=self.cluster_labels,
                plot_name="%sxyz.png" % self.params.plot_prefix,
            )
            _plot(
                (coord_reduced_x, coord_reduced_y, coord_reduced_z),
                labels=self.cluster_labels,
                plot_name="%sxyz_pca.png" % self.params.plot_prefix,
            )

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
            d["sym_op_scores"] = dict(
                (
                    (str(sym_op), score.as_dict())
                    for sym_op, score in self._symmetry_analysis.sym_op_scores.items()
                )
            )
            d["subgroup_scores"] = [
                score.as_dict() for score in self._symmetry_analysis.subgroup_scores
            ]
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
        import json

        json_str = json.dumps(d, indent=indent)
        if filename is not None:
            with open(filename, "wb") as f:
                f.write(json_str)
        return json.dumps(d, indent=indent)


def _plot(coords, labels=None, plot_name="xy.png"):
    from matplotlib import pyplot as plt
    import numpy

    coord_x = coords[0]
    coord_y = coords[1]

    fig = plt.figure()
    if len(coords) > 2:
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

        ax = fig.add_subplot(111, projection="3d")
        coord_z = coords[2]
    else:
        ax = fig.add_subplot(111)
        coord_z = None

    if labels is None:
        labels = flex.int(len(coord_x), -1)

    unique_labels = set(labels)
    unique_labels = sorted(unique_labels)
    n_clusters = len(unique_labels) - (1 if -1 in unique_labels else 0)
    colours = list(plt.cm.Spectral(numpy.linspace(0, 1, n_clusters)))
    if -1 in unique_labels:
        colours.insert(0, (255, 255, 255, 1))
    for k, col in zip(unique_labels, colours):
        isel = (labels == k).iselection()
        if k == -1:  # or len(class_members) < min_cluster_size:
            # Black used for noise.
            col = "k"
            col = "0.25"  # mid-grey
            markersize = 1
            marker = "+"
            alpha = 0.1
        else:
            markersize = 2
            marker = "o"
            alpha = 0.5
        edgecolor = col
        if coord_z is None:
            ax.scatter(
                coord_x.select(isel),
                coord_y.select(isel),
                s=markersize,
                marker=marker,
                c=col,
                edgecolor=edgecolor,
                alpha=alpha,
            )
            if k >= 0:
                # plot cluster centroid
                ax.scatter(
                    flex.mean(coord_x.select(isel)),
                    flex.mean(coord_y.select(isel)),
                    s=markersize * 10,
                    marker=marker,
                    c=col,
                    edgecolor="black",
                )
        else:
            ax.scatter(
                coord_x.select(isel),
                coord_y.select(isel),
                coord_z.select(isel),
                s=markersize,
                marker=marker,
                c=col,
                edgecolor=edgecolor,
                alpha=alpha,
            )
            if k >= 0:
                # plot cluster centroid
                ax.scatter(
                    flex.mean(coord_x.select(isel)),
                    flex.mean(coord_y.select(isel)),
                    flex.mean(coord_z.select(isel)),
                    s=markersize * 10,
                    marker=marker,
                    c=col,
                    edgecolor="black",
                )

    lim = max([1, flex.max(flex.abs(coord_x)), flex.max(flex.abs(coord_y))])
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    plt.savefig(plot_name, size_inches=(10, 10), dpi=300, bbox_inches="tight")
    plt.close(fig)


def _plot_angles(coords, labels=None, plot_name="phi_r.png"):
    coord_x, coord_y = coords

    r = flex.sqrt(flex.pow2(coord_x) + flex.pow2(coord_y))
    phi = flex.atan2(coord_y, coord_x)

    import math

    phi_deg = (180 / math.pi) * phi

    from matplotlib import pyplot as plt
    import numpy

    fig = plt.figure()
    ax = fig.add_subplot(111)

    if labels is None:
        labels = flex.int(len(coord_x), -1)

    unique_labels = set(labels)
    unique_labels = sorted(unique_labels)
    n_clusters = len(unique_labels) - (1 if -1 in unique_labels else 0)
    colours = list(plt.cm.Spectral(numpy.linspace(0, 1, n_clusters)))
    if -1 in unique_labels:
        colours.insert(0, (255, 255, 255, 1))
    for k, col in zip(unique_labels, colours):
        isel = (labels == k).iselection()
        if k == -1:  # or len(class_members) < min_cluster_size:
            # Black used for noise.
            col = "k"
            col = "0.25"  # mid-grey
            markersize = 1
            marker = "+"
        else:
            markersize = 2
            marker = "o"
        if 0 and not isinstance(col, basestring) and len(col) == 4:
            # darken the edges
            frac = 0.75
            edgecolor = [col[0] * frac, col[1] * frac, col[2] * frac, col[3]]
        else:
            edgecolor = col
        ax.scatter(
            phi_deg.select(isel),
            r.select(isel),
            s=markersize,
            marker=marker,
            c=col,
            edgecolor=edgecolor,
        )

    ax.set_xlim(-180, 180)
    ax.set_ylim(0, ax.get_ylim()[1])
    ax.set_xlabel("Angle ($^{\circ}$)")
    ax.set_ylabel("Magnitude")
    plt.savefig(plot_name, size_inches=(10, 10), dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_matrix(
    correlation_matrix, linkage_matrix, file_name, labels=None, color_threshold=0.05
):
    """Plot correlation and linkage matrices.

    Args:
      correlation_matrix (numpy.ndarray): The distance matrix used to generate
        the ``linkage_matrix``.
      linkage_matrix (numpy.ndarray): The hierarchical clustering of centroids of
        the initial clustering as produced by
        :func:`scipy.cluster.hierarchy.linkage`.
      file_name (str): The output file name.
      labels (list): Optional labels for the leaves of the dendrogram.
      color_threshold (float): The color threshold passed to the
        :func:scipy.cluster.hierarchy.dendrogram` function.

    """
    if correlation_matrix.shape[0] > 2000:
        return
    from matplotlib import pyplot as plt
    from scipy.cluster import hierarchy

    # Compute and plot dendrogram.
    fig = plt.figure(dpi=200, figsize=(16, 12))
    axdendro = fig.add_axes([0.09, 0.1, 0.2, 0.8])
    Y = linkage_matrix
    Z = hierarchy.dendrogram(Y, color_threshold=color_threshold, orientation="right")
    axdendro.set_xticks([])
    axdendro.set_yticks([])

    # Plot distance matrix.
    axmatrix = fig.add_axes([0.3, 0.1, 0.6, 0.8])
    index = Z["leaves"]
    D = correlation_matrix
    D = D[index, :]
    D = D[:, index]
    im = axmatrix.matshow(D, aspect="auto", origin="lower")
    axmatrix.yaxis.tick_right()
    if labels is not None:
        axmatrix.xaxis.tick_bottom()
        axmatrix.set_xticks(list(range(len(labels))))
        axmatrix.set_xticklabels([labels[i] for i in index], rotation=70)
        axmatrix.yaxis.set_ticks([])

    # Plot colorbar.
    axcolor = fig.add_axes([0.91, 0.1, 0.02, 0.8])
    plt.colorbar(im, cax=axcolor)

    # Display and save figure.
    fig.savefig(file_name)
    plt.close(fig)


def plot_dendrogram(linkage_matrix, file_name, labels=None, color_threshold=0.05):
    """Plot dendrogram from a linkage matrix.

    Args:
      linkage_matrix (numpy.ndarray): The hierarchical clustering of centroids of
        the initial clustering as produced by
        :func:`scipy.cluster.hierarchy.linkage`.
      file_name (str): The output file name.
      labels (list): Optional labels for the leaves of the dendrogram.
      color_threshold (float): The color threshold passed to the
        :func:scipy.cluster.hierarchy.dendrogram` function.

    """
    from matplotlib import pyplot as plt

    fig = plt.figure(dpi=200, figsize=(16, 12))

    from scipy.cluster import hierarchy

    ddict = hierarchy.dendrogram(
        linkage_matrix,
        color_threshold=color_threshold,
        labels=labels,
        show_leaf_counts=False,
    )
    locs, labels = plt.xticks()
    plt.setp(labels, rotation=70)
    fig.savefig(file_name)
    plt.close(fig)


class SymmetryAnalysis(object):
    def __init__(self, coords, sym_ops, subgroups, cb_op_inp_min):

        import scipy.spatial.distance as ssd

        self.subgroups = subgroups
        self.cb_op_inp_min = cb_op_inp_min
        lattice_group = subgroups.result_groups[0]["subsym"].space_group()
        X = coords.as_numpy_array()
        n_datasets = coords.all()[0] // len(sym_ops)
        dist_mat = ssd.pdist(X, metric="cosine")
        cos_angle = 1 - ssd.squareform(dist_mat)

        self._sym_ops_cos_angle = OrderedDict()
        for dataset_id in range(n_datasets):
            ref_sym_op_id = None
            ref_cluster_id = None
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
            ScoreSubGroup(subgrp, self.sym_op_scores.values())
            for subgrp in self.subgroups.result_groups
        ]
        total_likelihood = sum(score.likelihood for score in subgroup_scores)
        sort_order = flex.sort_permutation(
            flex.double(score.likelihood for score in subgroup_scores),
            reverse=True,
            stable=True,
        )
        self.subgroup_scores = [subgroup_scores[i] for i in sort_order]
        for score in self.subgroup_scores:
            score.likelihood /= total_likelihood

        # The 'confidence' scores are derived from the total probability of the best
        # solution p_best and that for the next best solution p_next:
        #   confidence = [p_best * (p_best - p_next)]^1/2.

        confidence = flex.double(len(self.subgroup_scores), 0)
        for i, score in enumerate(self.subgroup_scores[:-1]):
            next_score = self.subgroup_scores[i + 1]
            if score.likelihood > 0 and next_score.likelihood > 0:
                lgc = score.likelihood * (score.likelihood - next_score.likelihood)
                confidence = abs(lgc) ** 0.5
                if lgc < 0:
                    confidence = -confidence
                score.confidence = confidence

        self.best_solution = self.subgroup_scores[0]

    def __str__(self):
        """Return a string representation of the results.

        Returns:
          str:

        """
        output = []
        header = ("likelihood", "Z-CC", "CC", "", "Operator")
        rows = [header]
        for score in self.sym_op_scores.values():
            if score.likelihood > 0.9:
                stars = "***"
            elif score.likelihood > 0.7:
                stars = "**"
            elif score.likelihood > 0.5:
                stars = "*"
            else:
                stars = ""
            rows.append(
                (
                    "%.3f" % score.likelihood,
                    "%.2f" % score.z_cc,
                    "%.2f" % score.cc,
                    stars,
                    "%s" % score.sym_op.r().info(),
                )
            )
        output.append("Scoring individual symmetry elements")
        output.append(table_utils.format(rows, has_header=True, delim="  "))

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
        for score in self.subgroup_scores:
            if score.likelihood > 0.8:
                stars = "***"
            elif score.likelihood > 0.6:
                stars = "**"
            elif score.likelihood > 0.4:
                stars = "*"
            else:
                stars = ""
            rows.append(
                (
                    "%s" % score.subgroup["best_subsym"].space_group_info(),
                    stars,
                    "%.3f" % score.likelihood,
                    "% .2f" % score.z_cc_net,
                    "% .2f" % score.z_cc_for,
                    "% .2f" % score.z_cc_against,
                    "%.1f" % score.subgroup["max_angular_difference"],
                    "%s" % (score.subgroup["cb_op_inp_best"] * self.cb_op_inp_min),
                )
            )
        output.append("Scoring all possible sub-groups")
        output.append(table_utils.format(rows, has_header=True, delim="  "))

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
            "cc_unrelated_pairs": self.corr_unrelated.coefficient(),
            "n_unrelated_pairs": self.corr_unrelated.n(),
            "E_cc_true": self.E_cc_true,
            "cc_sig_fac": self.cc_sig_fac,
            "cc_true": self.cc_true,
        }

        d["sym_op_scores"] = [score.as_dict() for score in self.sym_op_scores]
        d["subgroup_scores"] = [score.as_dict() for score in self.subgroup_scores]
        return d


from dials.algorithms.symmetry.determine_space_group import ScoreCorrelationCoefficient


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
        return {"likelihood": self.likelihood, "z_cc": self.z_cc, "cc": self.cc}


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
            "likelihood": self.likelihood,
            "confidence": self.confidence,
            "z_cc_net": "% .2f" % self.z_cc_net,
            "z_cc_for": "% .2f" % self.z_cc_for,
            "z_cc_against": "% .2f" % self.z_cc_against,
            # "cc_for": "% .2f" % self.cc_for.coefficient(),
            # "cc_against": "% .2f" % self.cc_against.coefficient(),
            "max_angular_difference": "%.1f" % self.subgroup["max_angular_difference"],
            "cb_op": "%s" % (self.subgroup["cb_op_inp_best"]),
        }
