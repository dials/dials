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
from cctbx import sgtbx
import iotbx.phil

from dials.algorithms.symmetry.cosym import target
from dials.algorithms.symmetry.cosym import engine
from dials.algorithms.symmetry import symmetry_base

phil_scope = iotbx.phil.parse('''\

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
  agglomerative {
    n_clusters = 2
      .type = int(value_min=1)
  }
  seed {
    min_silhouette_score = 0.2
      .type = float(value_min=-1, value_max=1)
    n_clusters = auto
      .type = int(value_min=1)
  }
}

nproc = 1
  .type = int(value_min=1)
  .help = "The number of processes to use."

''')

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
      absolute_angle_tolerance=None)

    self.params = params
    self.intensities = self.intensities.customized_copy(
      space_group_info=self.input_space_group.change_basis(
        self.cb_op_inp_min).info())
    if self.params.dimensions is Auto:
      dimensions = None
    else:
      dimensions = self.params.dimensions
    lattice_group = None
    if self.params.lattice_group is not None:
      lattice_group = self.params.lattice_group.group() \
        .build_derived_patterson_group().info().primitive_setting().group()
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
        logger.info('Functional: %g' % self.minimizer.f)
        self._principal_component_analysis()
        dimensions.append(dim)
        functional.append(self.minimizer.f)
        explained_variance.append(self.explained_variance)
        explained_variance_ratio.append(self.explained_variance_ratio)

      # Find the elbow point of the curve, in the same manner as that used by
      # distl spotfinder for resolution method 1 (Zhang et al 2006).
      # See also dials/algorithms/spot_finding/per_image_analysis.py

      from scitbx import matrix
      x = flex.double(dimensions)
      y = flex.double(functional)
      slopes = (y[-1] - y[:-1])/(x[-1] - x[:-1])
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

      logger.info('Best number of dimensions: %i' %x_g)
      self.target.set_dimensions(int(x_g))

      if params.save_plot:
        from matplotlib import pyplot as plt
        fig = plt.figure(figsize=(10,8))
        plt.clf()
        plt.plot(dimensions, functional)
        plt.plot([x_g, x_g], plt.ylim())
        plt.xlabel('Dimensions')
        plt.ylabel('Functional')
        plt.savefig('%sfunctional_vs_dimension.png' % params.plot_prefix)

        plt.clf()
        for dim, expl_var in zip(dimensions, explained_variance):
          plt.plot(range(1, dim+1), expl_var, label='%s' % dim)
        plt.plot([x_g, x_g], plt.ylim())
        plt.xlabel('Dimension')
        plt.ylabel('Explained variance')
        plt.savefig('%sexplained_variance_vs_dimension.png' % params.plot_prefix)

        plt.clf()
        for dim, expl_var_ratio in zip(dimensions, explained_variance_ratio):
          plt.plot(range(1, dim+1), expl_var_ratio, label='%s' % dim)
        plt.plot([x_g, x_g], plt.ylim())
        plt.xlabel('Dimension')
        plt.ylabel('Explained variance ratio')
        plt.savefig(
          '%sexplained_variance_ratio_vs_dimension.png' % params.plot_prefix)
        plt.close(fig)

    self._optimise()
    self._principal_component_analysis()

    self._cosine_analysis()
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
      #min_iterations=tp.min_iterations,
      max_iterations=tp.max_iterations,
      max_calls=tp.max_calls
    )

    M = engine.lbfgs_with_curvs(
      self.target, coords,
      use_curvatures=self.params.use_curvatures,
      termination_params=termination_params
    )
    self.minimizer = M

    coords = M.x.deep_copy()
    coords.reshape(flex.grid(dim, NN*n_sym_ops))
    coords.matrix_transpose_in_place()
    self.coords = coords

  def _principal_component_analysis(self):
    # Perform PCA
    from sklearn.decomposition import PCA
    X = self.coords.as_numpy_array()
    pca = PCA().fit(X)
    logger.info('Principal component analysis:')
    logger.info('Explained variance: ' + ', '.join(
      ['%.2g' %v for v in pca.explained_variance_]))
    logger.info('Explained variance ratio: ' + ', '.join(
      ['%.2g' %v for v in pca.explained_variance_ratio_]))
    self.explained_variance = pca.explained_variance_
    self.explained_variance_ratio = pca.explained_variance_ratio_
    if self.target.dim > 3:
      pca.n_components = 3
    x_reduced = pca.fit_transform(X)

    import numpy
    self.coords_reduced = flex.double(numpy.ascontiguousarray(x_reduced))

  def _cosine_analysis(self):
    from scipy.cluster import hierarchy
    import scipy.spatial.distance as ssd

    X = self.coords.as_numpy_array()
    dist_mat = ssd.pdist(X, metric='cosine')
    cos_angle = 1 - ssd.squareform(dist_mat)
    linkage_matrix = hierarchy.linkage(dist_mat, method='average')

    c, coph_dists = hierarchy.cophenet(linkage_matrix, dist_mat)
    logger.debug(
      'Cophenetic correlation coefficient between heirarchical clustering and pairwise distance matrix: %.3f' % c)

    if self.params.save_plot:
      plot_matrix(
        cos_angle,
        linkage_matrix, '%scos_angle_matrix.png' % self.params.plot_prefix)
      plot_dendrogram(
        linkage_matrix,
        '%scos_angle_dendrogram.png' % self.params.plot_prefix)

    sym_ops = [sgtbx.rt_mx(s).new_denominators(1, 12) for s in self.target.get_sym_ops()]

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
    rows = [['symop', 'order', 'sg', 'mean(cos(angle))', 'median(cos(angle))']]
    perm = flex.sort_permutation(
      flex.double([flex.mean(ca) for ca in sym_ops_cos_angle.values()]),
      reverse=True)
    for p in perm:
      op, ca = sym_ops_cos_angle.items()[p]
      sg.expand_smx(op)
      rows.append((
        str(op), str(op.r().order()),
        str(sg.info().reference_setting()),
        '%.3f' % flex.mean(ca),
        '%.3f' % flex.median(ca),
      ))
    logger.info(
      'Analysis of cos(angle) between points corresponding to the same datasets:')
    logger.info(table_utils.format(rows, has_header=True))

  def _space_group_for_dataset(self, dataset_id, sym_ops):
    sg = copy.deepcopy(self.input_space_group)
    ref_sym_op_id = None
    ref_cluster_id = None
    for sym_op_id in range(len(sym_ops)):
      i_cluster = self.cluster_labels[len(self.input_intensities) * sym_op_id + dataset_id]
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
    n_clusters = len(set(self.cluster_labels)) - (1 if -1 in self.cluster_labels else 0)

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
              reindexing_ops[i_cluster] = partition[0].as_xyz()

    return reindexing_ops

  def _cluster_analysis(self):

    self.cluster_labels = self._do_clustering(self.params.cluster.method)

    cluster_miller_arrays = []

    space_groups = []

    sym_ops = [sgtbx.rt_mx(s).new_denominators(1, 12) for s in self.target.get_sym_ops()]
    self.space_groups = space_groups

    reindexing_ops = {}
    space_groups = {}

    for dataset_id in range(len(self.input_intensities)):
      space_groups[dataset_id] = self._space_group_for_dataset(
        dataset_id, sym_ops)

      cosets = sgtbx.cosets.left_decomposition(
        self.target._lattice_group,
        space_groups[dataset_id].info().primitive_setting().group())

      reindexing_ops[dataset_id] = self._reindexing_ops_for_dataset(
        dataset_id, sym_ops, cosets)

    self.space_groups = space_groups
    self.reindexing_ops = reindexing_ops

  def _do_clustering(self, method):
    if method == 'dbscan':
      clustering = self._dbscan_clustering
    elif method == 'bisect':
      clustering = self._bisect_clustering
    elif method == 'minimize_divide':
      clustering = self._minimize_divide_clustering
    elif method == 'agglomerative':
      clustering = self._agglomerative_clustering
    elif method == 'seed':
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
      min_samples=self.params.cluster.dbscan.min_samples
    ).fit(X)
    import numpy as np
    return flex.int(db.labels_.astype(np.int32))

  def _bisect_clustering(self):
    axis = self.params.cluster.bisect.axis
    assert axis < self.coords_reduced.all()[1]
    x = self.coords_reduced[:,axis:axis+1].as_1d()
    cluster_labels = flex.int(x.size(), 0)
    cluster_labels.set_selected(x > 0, 1)
    return cluster_labels

  def _minimize_divide_clustering(self):
    x = self.coords_reduced[:,:1].as_1d()
    y = self.coords_reduced[:,1:2].as_1d()
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
      n_clusters=self.params.cluster.agglomerative.n_clusters,
      linkage='average', affinity='cosine')
    model.fit(X)
    return flex.int(model.labels_.astype(np.int32))

  def _seed_clustering(self):
    from dials.algorithms.symmetry.cosym.seed_clustering import seed_clustering
    clustering = seed_clustering(
      self.coords,
      len(self.input_intensities),
      len(self.target.get_sym_ops()),
      min_silhouette_score=self.params.cluster.seed.min_silhouette_score,
      n_clusters=self.params.cluster.seed.n_clusters,
      plot_prefix=self.params.plot_prefix if self.params.save_plot else None
    )
    return clustering.cluster_labels

  def _plot(self):
    self.target.plot_rij_matrix(plot_name='%srij.png' % self.params.plot_prefix)
    self.target.plot_rij_histogram(
      plot_name='%srij_hist.png' % self.params.plot_prefix)
    self.target.plot_rij_cumulative_frequency(
      plot_name='%srij_sorted.png' % self.params.plot_prefix)
    self.target.plot_wij_matrix(plot_name='%swij.png' % self.params.plot_prefix)
    self.target.plot_wij_histogram(
      plot_name='%swij_hist.png' % self.params.plot_prefix)
    self.target.plot_wij_cumulative_frequency(
      plot_name='%swij_sorted.png' % self.params.plot_prefix)

    coord_x = self.coords[:,0:1].as_1d()
    coord_y = self.coords[:,1:2].as_1d()
    assert coord_x.size() == coord_y.size(), (coord_x.size(), coord_y.size())
    coord_reduced_x = self.coords_reduced[:,0:1].as_1d()
    coord_reduced_y = self.coords_reduced[:,1:2].as_1d()
    _plot((coord_x, coord_y), labels=self.cluster_labels,
          plot_name='%sxy.png' % self.params.plot_prefix)
    _plot((coord_reduced_x, coord_reduced_y), labels=self.cluster_labels,
          plot_name='%sxy_pca.png' % self.params.plot_prefix)
    _plot_angles((coord_x, coord_y), labels=self.cluster_labels,
                 plot_name='%sphi_r.png' % self.params.plot_prefix)
    _plot_angles((coord_reduced_x, coord_reduced_y), labels=self.cluster_labels,
                 plot_name='%sphi_r_pca.png' % self.params.plot_prefix)

    if self.coords_reduced.all()[1] > 2:
      coord_z = self.coords[:,2:3].as_1d()
      coord_reduced_z = self.coords_reduced[:,2:3].as_1d()
      _plot((coord_x, coord_y, coord_z), labels=self.cluster_labels,
            plot_name='%sxyz.png' % self.params.plot_prefix)
      _plot((coord_reduced_x, coord_reduced_y, coord_reduced_z),
            labels=self.cluster_labels,
            plot_name='%sxyz_pca.png' % self.params.plot_prefix)

def _plot(coords, labels=None, plot_name='xy.png'):
  from matplotlib import pyplot as plt
  import numpy

  coord_x = coords[0]
  coord_y = coords[1]

  fig = plt.figure()
  if len(coords) > 2:
    from mpl_toolkits.mplot3d import Axes3D # import dependency
    ax = fig.add_subplot(111, projection='3d')
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
    colours.insert(0, (255,255,255,1))
  for k, col in zip(unique_labels, colours):
    isel = (labels == k).iselection()
    if k == -1:# or len(class_members) < min_cluster_size:
        # Black used for noise.
      col = 'k'
      col = '0.25' # mid-grey
      markersize = 1
      marker = '+'
      alpha = 0.1
    else:
      markersize = 2
      marker = 'o'
      alpha = 0.5
    edgecolor = col
    if coord_z is None:
      ax.scatter(coord_x.select(isel), coord_y.select(isel),
                 s=markersize, marker=marker, c=col, edgecolor=edgecolor, alpha=alpha)
      if k >= 0:
        # plot cluster centroid
        ax.scatter(flex.mean(coord_x.select(isel)), flex.mean(coord_y.select(isel)),
                   s=markersize*10, marker=marker, c=col, edgecolor='black')
    else:
      ax.scatter(coord_x.select(isel), coord_y.select(isel), coord_z.select(isel),
                 s=markersize, marker=marker, c=col, edgecolor=edgecolor, alpha=alpha)
      if k >= 0:
        # plot cluster centroid
        ax.scatter(flex.mean(coord_x.select(isel)), flex.mean(coord_y.select(isel)),
                   flex.mean(coord_z.select(isel)),
                   s=markersize*10, marker=marker, c=col, edgecolor='black')

  ax.set_xlim(-1,1)
  ax.set_ylim(-1,1)
  ax.set_aspect("equal")
  plt.savefig(plot_name,
              size_inches=(10,10),
              dpi=300,
              bbox_inches='tight')
  plt.close(fig)


def _plot_angles(coords, labels=None, plot_name='phi_r.png'):
  coord_x, coord_y = coords

  r = flex.sqrt(flex.pow2(coord_x) + flex.pow2(coord_y))
  phi = flex.atan2(coord_y, coord_x)

  import math
  phi_deg = (180/math.pi) * phi

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
    colours.insert(0, (255,255,255,1))
  for k, col in zip(unique_labels, colours):
    isel = (labels == k).iselection()
    if k == -1:# or len(class_members) < min_cluster_size:
        # Black used for noise.
      col = 'k'
      col = '0.25' # mid-grey
      markersize = 1
      marker = '+'
    else:
      markersize = 2
      marker = 'o'
    if 0 and not isinstance(col, basestring) and len(col) == 4:
      # darken the edges
      frac = 0.75
      edgecolor = [col[0]*frac, col[1]*frac, col[2]*frac, col[3]]
    else:
      edgecolor = col
    ax.scatter(phi_deg.select(isel), r.select(isel),
               s=markersize, marker=marker, c=col, edgecolor=edgecolor)

  ax.set_xlim(-180,180)
  ax.set_ylim(0,ax.get_ylim()[1])
  ax.set_xlabel('Angle ($^{\circ}$)')
  ax.set_ylabel('Magnitude')
  plt.savefig(plot_name,
              size_inches=(10,10),
              dpi=300,
              bbox_inches='tight')
  plt.close(fig)


def plot_matrix(correlation_matrix, linkage_matrix, file_name, labels=None,
                color_threshold=0.05):
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
  fig = plt.figure(dpi=200, figsize=(16,12))
  axdendro = fig.add_axes([0.09,0.1,0.2,0.8])
  Y = linkage_matrix
  Z = hierarchy.dendrogram(Y,
                           color_threshold=color_threshold,
                           orientation='right')
  axdendro.set_xticks([])
  axdendro.set_yticks([])

  # Plot distance matrix.
  axmatrix = fig.add_axes([0.3,0.1,0.6,0.8])
  index = Z['leaves']
  D = correlation_matrix
  D = D[index,:]
  D = D[:,index]
  im = axmatrix.matshow(D, aspect='auto', origin='lower')
  axmatrix.yaxis.tick_right()
  if labels is not None:
    axmatrix.xaxis.tick_bottom()
    axmatrix.set_xticks(list(range(len(labels))))
    axmatrix.set_xticklabels([labels[i] for i in index], rotation=70)
    axmatrix.yaxis.set_ticks([])

  # Plot colorbar.
  axcolor = fig.add_axes([0.91,0.1,0.02,0.8])
  plt.colorbar(im, cax=axcolor)

  # Display and save figure.
  fig.savefig(file_name)
  plt.close(fig)


def plot_dendrogram(linkage_matrix, file_name, labels=None,
                    color_threshold=0.05):
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

  fig = plt.figure(dpi=200, figsize=(16,12))

  from scipy.cluster import hierarchy
  ddict = hierarchy.dendrogram(linkage_matrix,
                               color_threshold=color_threshold,
                               labels=labels,
                               show_leaf_counts=False)
  locs, labels = plt.xticks()
  plt.setp(labels, rotation=70)
  fig.savefig(file_name)
  plt.close(fig)
