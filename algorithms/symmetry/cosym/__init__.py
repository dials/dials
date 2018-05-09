from __future__ import absolute_import, division, print_function

import logging
logger = logging.getLogger(__name__)

import copy
from collections import OrderedDict
import math

from libtbx import Auto
from libtbx.utils import Sorry
from libtbx import table_utils
from scitbx.array_family import flex
from cctbx import sgtbx
import iotbx.phil

from dials.algorithms.symmetry.cosym import target
from dials.algorithms.symmetry.cosym import engine

phil_scope = iotbx.phil.parse('''\

lattice_group = None
  .type = space_group

animate = False
  .type = bool

save_intermediate_plots = False
  .type = bool

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

verbose = False
  .type = bool

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

class analyse_datasets(object):

  def __init__(self, datasets, params):
    self.datasets = datasets
    self.params = params

    self.input_space_group = None
    for dataset in datasets:
      if self.input_space_group is None:
        self.input_space_group = dataset.space_group()
      else:
        assert dataset.space_group() == self.input_space_group

    if self.params.dimensions is Auto:
      dimensions = None
    else:
      dimensions = self.params.dimensions
    lattice_group = None
    if self.params.lattice_group is not None:
      lattice_group = self.params.lattice_group.group()
    self.target = target.Target(
      self.datasets,
      min_pairs=self.params.min_pairs,
      lattice_group=lattice_group,
      dimensions=dimensions,
      verbose=self.params.verbose,
      weights=self.params.weights,
      nproc=self.params.nproc,
    )
    if self.params.dimensions is Auto:
      dimensions = []
      functional = []
      explained_variance = []
      explained_variance_ratio = []
      for dim in range(1, self.target.dim + 1):
        self.target.set_dimensions(dim)
        self.optimise()
        logger.info('Functional: %g' % self.minimizer.f)
        self.principal_component_analysis()
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

    self.optimise()
    self.principal_component_analysis()

    self.cosine_analysis()
    self.cluster_analysis()
    if self.params.save_plot:
      self.plot()

  def optimise(self):

    NN = len(self.datasets)
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
      self.target, coords, verbose=True, animate=self.params.animate,
      save_intermediate_plots=self.params.save_intermediate_plots,
      use_curvatures=self.params.use_curvatures,
      termination_params=termination_params
    )
    self.minimizer = M

    coords = M.x.deep_copy()
    coords.reshape(flex.grid(dim, NN*n_sym_ops))
    coords.matrix_transpose_in_place()
    self.coords = coords

  def principal_component_analysis(self):
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

  def cosine_analysis(self):
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

    for dataset_id in range(len(self.datasets)):
      ref_sym_op_id = None
      ref_cluster_id = None
      for sym_op_id in range(len(sym_ops)):
        if ref_sym_op_id is None:
          ref_sym_op_id = sym_op_id
          continue
        op = sym_ops[ref_sym_op_id].inverse().multiply(sym_ops[sym_op_id])
        op = op.new_denominators(1, 12)

        ref_idx = len(self.datasets) * ref_sym_op_id + dataset_id
        comp_idx = len(self.datasets) * sym_op_id + dataset_id
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

  def cluster_analysis(self):
    from cctbx.sgtbx import cosets

    if self.params.cluster.method == 'dbscan':
      self.dbscan_clustering()
    elif self.params.cluster.method == 'bisect':
      self.bisect_clustering()
    elif self.params.cluster.method == 'minimize_divide':
      self.minimize_divide_clustering()
    elif self.params.cluster.method == 'agglomerative':
      self.agglomerative_clustering()
    elif self.params.cluster.method == 'seed':
      self.seed_clustering()

    # Number of clusters in labels, ignoring noise if present.
    n_clusters = len(set(self.cluster_labels)) - (1 if -1 in self.cluster_labels else 0)

    cluster_miller_arrays = []

    space_groups = []

    sym_ops = [sgtbx.rt_mx(s).new_denominators(1, 12) for s in self.target.get_sym_ops()]
    self.space_groups = space_groups

    reindexing_ops = {}
    space_groups = {}

    for dataset_id in range(len(self.datasets)):
      sg = copy.deepcopy(self.input_space_group)
      ref_sym_op_id = None
      ref_cluster_id = None
      for sym_op_id in range(len(sym_ops)):
        i_cluster = self.cluster_labels[len(self.datasets) * sym_op_id + dataset_id]
        if i_cluster < 0:
          continue
        if ref_sym_op_id is None:
          ref_sym_op_id = sym_op_id
          ref_cluster_id = i_cluster
          continue
        op = sym_ops[ref_sym_op_id].inverse().multiply(sym_ops[sym_op_id])
        if i_cluster == ref_cluster_id:
          sg.expand_smx(op.new_denominators(1, 12))

      sg.make_tidy()
      space_groups[dataset_id] = sg

      coset = cosets.left_decomposition(
        self.target._lattice_group,
        sg.info().primitive_setting().group())

      reindexing_ops[dataset_id] = {}

      for i_cluster in range(n_clusters):
        isel = (self.cluster_labels == i_cluster).iselection()
        dataset_ids = isel % len(self.datasets)
        idx = flex.first_index(dataset_ids, dataset_id)
        sel = (dataset_ids == dataset_id).iselection()
        if idx >= 0:
          sym_op_id = isel[idx] // len(self.datasets)
        for s in sel:
          sym_op_id = isel[s] // len(self.datasets)

          for partition in coset.partitions:
            if sym_ops[sym_op_id] in partition:
              if i_cluster not in reindexing_ops[dataset_id]:
                reindexing_ops[dataset_id][i_cluster] = partition[0].as_xyz()
              #else:
                #assert reindexing_ops[dataset_id][i_cluster] == partition[0].as_xyz()

    self.space_groups = space_groups
    self.reindexing_ops = reindexing_ops

  def dbscan_clustering(self):
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
    self.cluster_labels = flex.int(db.labels_.astype(np.int32))

  def bisect_clustering(self):
    axis = self.params.cluster.bisect.axis
    assert axis < self.coords_reduced.all()[1]
    x = self.coords_reduced[:,axis:axis+1].as_1d()
    self.cluster_labels = flex.int(x.size(), 0)
    self.cluster_labels.set_selected(x > 0, 1)

  def minimize_divide_clustering(self):
    x = self.coords_reduced[:,:1].as_1d()
    y = self.coords_reduced[:,1:2].as_1d()
    from cctbx.merging.brehm_diederichs import minimize_divide
    selection = minimize_divide(x, y).plus_minus()
    self.cluster_labels = flex.int(x.size(), 0)
    self.cluster_labels.set_selected(selection, 1)

  def agglomerative_clustering(self):
    X = self.coords.as_numpy_array()

    # Perform cluster analysis
    from sklearn.cluster import AgglomerativeClustering
    import numpy as np
    model = AgglomerativeClustering(
      n_clusters=self.params.cluster.agglomerative.n_clusters,
      linkage='average', affinity='cosine')
    model.fit(X)
    self.cluster_labels = flex.int(model.labels_.astype(np.int32))

  def seed_clustering(self):
    eps = 1e-6
    X_orig = self.coords.as_numpy_array()

    import numpy as np
    from scipy.cluster import hierarchy
    import scipy.spatial.distance as ssd
    from sklearn.neighbors import NearestNeighbors
    from sklearn import metrics

    # initialise cluster labels: -1 signifies doesn't belong to a cluster
    self.cluster_labels = flex.int(self.coords.all()[0], -1)

    cluster_id = 0
    while self.cluster_labels.count(-1) > 0:
      dataset_ids = (flex.int_range(
        len(self.datasets) * len(self.target.get_sym_ops())) % len(self.datasets)
                     ).as_numpy_array()
      coord_ids = flex.int_range(dataset_ids.size).as_numpy_array()

      # select only those points that don't already belong to a cluster
      sel = np.where(self.cluster_labels == -1)
      X = X_orig[sel]
      dataset_ids = dataset_ids[sel]
      coord_ids = coord_ids[sel]

      # choose a high density point as seed for cluster
      nbrs = NearestNeighbors(
        n_neighbors=min(11, len(X)), algorithm='brute', metric='cosine').fit(X)
      distances, indices = nbrs.kneighbors(X)
      average_distance = flex.double([dist[1:].mean() for dist in distances])
      i = flex.min_index(average_distance)

      d_id = dataset_ids[i]
      cluster = np.array([coord_ids[i]])
      cluster_dataset_ids = np.array([d_id])
      xis = np.array([X[i]])

      for j in range(len(self.datasets)-1):
        # select only those rows that don't correspond to a dataset already
        # present in current cluster
        sel = np.where(dataset_ids != d_id)
        X = X[sel]
        dataset_ids = dataset_ids[sel]
        coord_ids = coord_ids[sel]

        assert len(X) > 0

        # Find nearest neighbour in cosine-space to the current cluster centroid
        nbrs = NearestNeighbors(
          n_neighbors=min(1, len(X)), algorithm='brute', metric='cosine').fit(X)
        distances, indices = nbrs.kneighbors([xis.mean(axis=0)])
        k = indices[0][0]
        d_id = dataset_ids[k]
        cluster = np.append(cluster, coord_ids[k])
        cluster_dataset_ids = np.append(cluster_dataset_ids, d_id)
        xis = np.append(xis, [X[k]], axis=0)

      # label this cluster
      self.cluster_labels.set_selected(flex.size_t(cluster.tolist()), cluster_id)
      cluster_id += 1

    if flex.max(self.cluster_labels) == 0:
      # assume single cluster
      return self.cluster_labels

    cluster_centroids = []
    X = self.coords.as_numpy_array()
    for i in set(self.cluster_labels):
      sel = self.cluster_labels == i
      cluster_centroids.append(
        X[(self.cluster_labels == i).iselection().as_numpy_array()].mean(axis=0))

    # hierarchical clustering of cluster centroids, using cosine metric
    dist_mat = ssd.pdist(cluster_centroids, metric='cosine')
    linkage_matrix = hierarchy.linkage(dist_mat, method='average')

    # compare valid equal-sized clustering using silhouette scores
    # https://en.wikipedia.org/wiki/Silhouette_(clustering)
    # http://scikit-learn.org/stable/auto_examples/cluster/plot_kmeans_silhouette_analysis.html
    distances = linkage_matrix[::, 2]
    distances = np.insert(distances, 0, 0)
    silhouette_scores = flex.double()
    thresholds = flex.double()
    n_clusters = flex.size_t()
    for threshold in distances[1:]:
      cluster_labels = self.cluster_labels.deep_copy()
      labels = hierarchy.fcluster(
        linkage_matrix, threshold-eps, criterion='distance').tolist()
      counts = [labels.count(l) for l in set(labels)]
      if len(set(counts)) > 1:
        # only equal-sized clusters are valid
        continue

      n = len(set(labels))
      if n == 1: continue
      for i in range(len(labels)):
        cluster_labels.set_selected(self.cluster_labels == i, int(labels[i]-1))
      silhouette_avg = metrics.silhouette_score(
        X, cluster_labels.as_numpy_array(), metric='cosine')
      # Compute the silhouette scores for each sample
      sample_silhouette_values = metrics.silhouette_samples(
        X, cluster_labels.as_numpy_array(), metric='cosine')
      silhouette_avg = sample_silhouette_values.mean()
      silhouette_scores.append(silhouette_avg)
      thresholds.append(threshold)
      n_clusters.append(n)

      count_negative = (sample_silhouette_values < 0).sum()
      logger.info('Clustering:')
      logger.info('  Number of clusters: %i' % n)
      logger.info('  Threshold score: %.3f (%.1f deg)' % (
        threshold, math.degrees(math.acos(1-threshold))))
      logger.info('  Silhouette score: %.3f' % silhouette_avg)
      logger.info('  -ve silhouette scores: %.1f%%' % (
        100 * count_negative/sample_silhouette_values.size))

      if self.params.save_plot:
        plot_silhouette(
          sample_silhouette_values, cluster_labels.as_numpy_array(),
          file_name='%ssilhouette_%i.png' % (self.params.plot_prefix, n))

    if self.params.cluster.seed.n_clusters is Auto:
      idx = flex.max_index(silhouette_scores)
    else:
      idx = flex.first_index(n_clusters, self.params.cluster.seed.n_clusters)
      if idx is None:
        raise Sorry('No valid clustering with %i clusters'
                    % self.params.cluster.seed.n_clusters)

    if (self.params.cluster.seed.n_clusters is Auto and
        silhouette_scores[idx] < self.params.cluster.seed.min_silhouette_score):
      # assume single cluster
      self.cluster_labels = flex.int(self.cluster_labels.size(), 0)
    else:
      threshold = thresholds[idx] - eps
      labels = hierarchy.fcluster(
        linkage_matrix, threshold, criterion='distance')
      cluster_labels = flex.double(self.cluster_labels.size(), -1)
      for i in range(len(labels)):
        cluster_labels.set_selected(self.cluster_labels == i, labels[i]-1)
      self.cluster_labels = cluster_labels

    if self.params.save_plot:
      plot_matrix(
        1 - ssd.squareform(dist_mat), linkage_matrix,
        '%sseed_clustering_cos_angle_matrix.png' % self.params.plot_prefix,
        color_threshold=threshold)
      plot_dendrogram(
        linkage_matrix,
        '%sseed_clustering_cos_angle_dendrogram.png' % self.params.plot_prefix,
        color_threshold=threshold)

    return self.cluster_labels

  def plot(self):
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
    coord_reduced_x = self.coords_reduced[:,0:1].as_1d()
    coord_reduced_y = self.coords_reduced[:,1:2].as_1d()
    plot((coord_x, coord_y), labels=self.cluster_labels,
         plot_name='%sxy.png' % self.params.plot_prefix)
    plot((coord_reduced_x, coord_reduced_y), labels=self.cluster_labels,
         plot_name='%sxy_pca.png' % self.params.plot_prefix)
    plot_angles((coord_x, coord_y), labels=self.cluster_labels,
                plot_name='%sphi_r.png' % self.params.plot_prefix)
    plot_angles((coord_reduced_x, coord_reduced_y), labels=self.cluster_labels,
                plot_name='%sphi_r_pca.png' % self.params.plot_prefix)

    if self.coords_reduced.all()[1] > 2:
      coord_z = self.coords[:,2:3].as_1d()
      coord_reduced_z = self.coords_reduced[:,2:3].as_1d()
      plot((coord_x, coord_y, coord_z), labels=self.cluster_labels,
           plot_name='%sxyz.png' % self.params.plot_prefix)
      plot((coord_reduced_x, coord_reduced_y, coord_reduced_z),
           labels=self.cluster_labels,
           plot_name='%sxyz_pca.png' % self.params.plot_prefix)

def plot(coords, labels=None, plot_centroids=True, plot_all=True, show=False, plot_name=None):
  assert plot_centroids or plot_all
  assert len(coords) >= 2

  coord_x = coords[0]
  coord_y = coords[1]
  if len(coords) > 2:
    coord_z = coords[2]
  else:
    coord_z = None

  from matplotlib import pyplot as plt
  import numpy

  fig = plt.figure()

  if coord_z is not None:
    from mpl_toolkits.mplot3d import Axes3D # import dependency
    ax = fig.add_subplot(111, projection='3d')
  else:
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
      alpha = 0.1
      #continue
    else:
      markersize = 2
      marker = 'o'
      alpha = 0.5
    if 0 and not isinstance(col, basestring) and len(col) == 4:
      # darken the edges
      frac = 0.75
      edgecolor = [col[0]*frac, col[1]*frac, col[2]*frac, col[3]]
    else:
      edgecolor = col
    if coord_z is None:
      if plot_all:
        ax.scatter(coord_x.select(isel), coord_y.select(isel),
                   s=markersize, marker=marker, c=col, edgecolor=edgecolor, alpha=alpha)
      if plot_centroids and k >= 0:
        ax.scatter(flex.mean(coord_x.select(isel)), flex.mean(coord_y.select(isel)),
                   s=markersize*10, marker=marker, c=col, edgecolor='black')
    else:
      if plot_all:
        ax.scatter(coord_x.select(isel), coord_y.select(isel), coord_z.select(isel),
                   s=markersize, marker=marker, c=col, edgecolor=edgecolor, alpha=alpha)
      if plot_centroids and k >= 0:
        ax.scatter(flex.mean(coord_x.select(isel)), flex.mean(coord_y.select(isel)),
                   flex.mean(coord_z.select(isel)),
                   s=markersize*10, marker=marker, c=col, edgecolor='black')

  ax.set_xlim(-1,1)
  ax.set_ylim(-1,1)
  ax.set_aspect("equal")
  if plot_name is not None:
    plt.savefig(plot_name,
                size_inches=(10,10),
                dpi=300,
                bbox_inches='tight')
  if show:
    plt.show()
  plt.close(fig)


def plot_angles(coords, labels=None, show=False, plot_name=None):
  assert len(coords) >= 2

  coord_x = coords[0]
  coord_y = coords[1]

  if len(coords) > 2:
    coord_z = coords[2]
  else:
    coord_z = None

  r = flex.sqrt(flex.pow2(coord_x) + flex.pow2(coord_y))
  phi = flex.atan2(coord_y, coord_x)

  import math
  phi_deg = (180/math.pi) * phi

  from matplotlib import pyplot as plt
  import numpy

  fig = plt.figure()

  if 0 and coord_z is not None:
    from mpl_toolkits.mplot3d import Axes3D # import dependency
    ax = fig.add_subplot(111, projection='3d')
  else:
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
      #continue
    else:
      markersize = 2
      marker = 'o'
    if 0 and not isinstance(col, basestring) and len(col) == 4:
      # darken the edges
      frac = 0.75
      edgecolor = [col[0]*frac, col[1]*frac, col[2]*frac, col[3]]
    else:
      edgecolor = col
    if coord_z is None:
      ax.scatter(phi_deg.select(isel), r.select(isel),
                 s=markersize, marker=marker, c=col, edgecolor=edgecolor)
    else:
      ax.scatter(phi_deg.select(isel), r.select(isel), coord_z.select(isel),
                 s=markersize, marker=marker, c=col, edgecolor=edgecolor)

  ax.set_xlim(-180,180)
  ax.set_ylim(0,ax.get_ylim()[1])
  ax.set_xlabel('Angle ($^{\circ}$)')
  ax.set_ylabel('Magnitude')
  if plot_name is not None:
    plt.savefig(plot_name,
                size_inches=(10,10),
                dpi=300,
                bbox_inches='tight')
  if show:
    plt.show()
  plt.close(fig)


def plot_matrix(correlation_matrix, linkage_matrix, file_name, labels=None,
                color_threshold=0.05):
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


def plot_silhouette(sample_silhouette_values, cluster_labels, file_name):
  import numpy as np
  from matplotlib import pyplot as plt
  import matplotlib.cm as cm

  fig = plt.figure()
  ax1 = fig.gca()
  y_lower = 10

  n_clusters = len(set(cluster_labels))
  silhouette_avg = sample_silhouette_values.mean()
  for i in range(n_clusters):
    # Aggregate the silhouette scores for samples belonging to
    # cluster i, and sort them
    ith_cluster_silhouette_values = \
      sample_silhouette_values[cluster_labels == i]

    ith_cluster_silhouette_values.sort()

    size_cluster_i = ith_cluster_silhouette_values.shape[0]
    y_upper = y_lower + size_cluster_i

    color = cm.spectral(float(i) / n_clusters)
    ax1.fill_betweenx(np.arange(y_lower, y_upper),
                      0, ith_cluster_silhouette_values,
                      facecolor=color, edgecolor=color, alpha=0.7)

    # Label the silhouette plots with their cluster numbers at the middle
    ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

    # Compute the new y_lower for next plot
    y_lower = y_upper + 10  # 10 for the 0 samples

  ax1.set_title("The silhouette plot for the various clusters.")
  ax1.set_xlabel("The silhouette coefficient values")
  ax1.set_ylabel("Cluster label")

  # The vertical line for average silhouette score of all the values
  ax1.axvline(x=silhouette_avg, color="red", linestyle="--")

  ax1.set_yticks([])  # Clear the yaxis labels / ticks
  ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])
  fig.savefig(file_name)
  plt.close(fig)
