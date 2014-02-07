#!/usr/bin/env python
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# dials.algorithms.indexing.fft3d.py
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Richard Gildea
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
import math

import libtbx
from scitbx import fftpack
from scitbx import matrix
from cctbx import crystal, uctbx, xray
from cctbx.array_family import flex
from dials.algorithms.indexing.indexer2 \
     import indexer_base, optimise_basis_vectors
from dials.algorithms.indexing.indexer2 \
     import vector_group, is_approximate_integer_multiple
from dials.model.experiment.experiment_list import Experiment, ExperimentList

class indexer_fft3d(indexer_base):

  def __init__(self, reflections, sweep, params):
    super(indexer_fft3d, self).__init__(reflections, sweep, params)

  def find_lattices(self):
    if self.params.reciprocal_space_grid.d_min is libtbx.Auto:
      # rough calculation of suitable d_min based on max cell
      # see also Campbell, J. (1998). J. Appl. Cryst., 31(3), 407-413.
      self.params.reciprocal_space_grid.d_min = (
        5 * self.params.max_cell / self.params.reciprocal_space_grid.n_points)
      print "Setting d_min: %s" %self.params.reciprocal_space_grid.d_min
    n_points = self.params.reciprocal_space_grid.n_points
    self.gridding = fftpack.adjust_gridding_triple(
      (n_points,n_points,n_points), max_prime=5)
    n_points = self.gridding[0]
    self.map_centroids_to_reciprocal_space_grid()
    self.d_min = self.params.reciprocal_space_grid.d_min

    print "Number of centroids used: %i" %(
      (self.reciprocal_space_grid>0).count(True))
    self.fft()
    if self.params.debug:
      self.debug_write_ccp4_map(map_data=self.grid_real, file_name="patt.map")
    self.find_peaks()
    if self.params.multiple_lattice_search:
      self.find_basis_vector_combinations_cluster_analysis()
      if self.params.debug:
        self.debug_show_candidate_basis_vectors()
      crystal_models = self.candidate_crystal_models
      if self.params.max_lattices is not None:
        crystal_models = crystal_models[:self.params.max_lattices]
    else:
      self.find_candidate_basis_vectors()
      # self.find_candidate_basis_vectors_nks()
      if self.params.debug:
        self.debug_show_candidate_basis_vectors()
      self.candidate_crystal_models = self.find_candidate_orientation_matrices(
        self.candidate_basis_vectors)
      crystal_models = self.candidate_crystal_models[:1]
    if self.target_symmetry_primitive is not None:
      crystal_models = [
        self.apply_symmetry(cm, self.target_symmetry_primitive)
        for cm in crystal_models]
    experiments = ExperimentList()
    for cm in crystal_models:
      experiments.append(Experiment(beam=self.beam,
                                    detector=self.detector,
                                    goniometer=self.goniometer,
                                    scan=self.scan,
                                    crystal=cm))
    return experiments

  def map_centroids_to_reciprocal_space_grid(self):
    self._map_to_grid_timer.start()
    assert len(self.reciprocal_space_points) == len(self.reflections)
    wavelength = self.beam.get_wavelength()
    d_min = self.params.reciprocal_space_grid.d_min

    n_points = self.gridding[0]
    rlgrid = 2 / (d_min * n_points)

    # real space FFT grid dimensions
    cell_lengths = [n_points * d_min/2 for i in range(3)]
    self.fft_cell = uctbx.unit_cell(cell_lengths+[90]*3)
    self.crystal_symmetry = crystal.symmetry(unit_cell=self.fft_cell,
                                             space_group_symbol="P1")

    print "FFT gridding: (%i,%i,%i)" %self.gridding

    grid = flex.double(flex.grid(self.gridding), 0)

    reflections_used_for_indexing = flex.size_t()

    for i_ref, point in enumerate(self.reciprocal_space_points):
      if self.reflections[i_ref].crystal != -1:
        continue
      point = matrix.col(point)
      spot_resolution = 1/point.length()
      if spot_resolution < d_min:
        continue

      grid_coordinates = [int(round(point[i]/rlgrid)+n_points/2) for i in range(3)]
      if max(grid_coordinates) >= n_points: continue # this reflection is outside the grid
      if min(grid_coordinates) < 0: continue # this reflection is outside the grid
      T = math.exp(-self.params.b_iso * point.length()**2 / 4)
      grid[grid_coordinates] = T
      reflections_used_for_indexing.append(i_ref)

    self.reciprocal_space_grid = grid
    self.reflections_used_for_indexing = reflections_used_for_indexing
    self._map_to_grid_timer.stop()

  def fft(self):
    self._fft_timer.start()
    #gb_to_bytes = 1073741824
    #bytes_to_gb = 1/gb_to_bytes
    #(128**3)*8*2*bytes_to_gb
    #0.03125
    #(256**3)*8*2*bytes_to_gb
    #0.25
    #(512**3)*8*2*bytes_to_gb
    #2.0

    fft = fftpack.complex_to_complex_3d(self.gridding)
    grid_complex = flex.complex_double(
      reals=self.reciprocal_space_grid,
      imags=flex.double(self.reciprocal_space_grid.size(), 0))
    grid_transformed = fft.forward(grid_complex)
    #self.grid_real = flex.pow2(flex.abs(grid_transformed))
    self.grid_real = flex.pow2(flex.real(grid_transformed))
    #self.grid_real = flex.pow2(flex.imag(self.grid_transformed))
    del grid_transformed
    self._fft_timer.stop()

  def find_peaks(self):
    self._find_peaks_timer.start()
    grid_real_binary = self.grid_real.deep_copy()
    rmsd = math.sqrt(
      flex.mean(flex.pow2(grid_real_binary.as_1d()-flex.mean(grid_real_binary.as_1d()))))
    grid_real_binary.set_selected(grid_real_binary < (self.params.rmsd_cutoff)*rmsd, 0)
    grid_real_binary.as_1d().set_selected(grid_real_binary.as_1d() > 0, 1)
    grid_real_binary = grid_real_binary.iround()
    from cctbx import masks
    flood_fill = masks.flood_fill(grid_real_binary, self.fft_cell)
    # the peak at the origin might have a significantly larger volume than the
    # rest so exclude this peak from determining maximum volume
    isel = (flood_fill.grid_points_per_void() > int(
        0.2 * flex.max(flood_fill.grid_points_per_void()[1:]))).iselection()
    #sites = flex.double()
    #for i in isel:
      #sel = grid_real_binary == (i+2)
      #max_value = flex.max(grid_real.as_1d().select(sel.as_1d()))
      #idx = flex.first_index(grid_real, max_value)
      #sites.append(

    if self.params.optimise_initial_basis_vectors:
      sites_cart = flood_fill.centres_of_mass_cart().select(isel)
      sites_cart_optimised = optimise_basis_vectors(
        #self.reciprocal_space_points,
        self.reciprocal_space_points.select(self.reflections_used_for_indexing),
        sites_cart)

      self.sites = self.fft_cell.fractionalize(sites_cart_optimised)

      diffs = (sites_cart_optimised - sites_cart)
      norms = diffs.norms()
      flex.min_max_mean_double(norms).show()
      perm = flex.sort_permutation(norms, reverse=True)
      for p in perm[:10]:
        print sites_cart[p], sites_cart_optimised[p], norms[p]

      # only use those vectors which haven't shifted too far from starting point
      self.sites = self.sites.select(
        norms < (5 * self.fft_cell.parameters()[0]/self.gridding[0]))
      #diff = (self.sites - flood_fill.centres_of_mass_frac().select(isel))
      #flex.min_max_mean_double(diff.norms()).show()

    else:
      self.sites = flood_fill.centres_of_mass_frac().select(isel)

    self._find_peaks_timer.stop()

  def find_candidate_basis_vectors(self):
    # hijack the xray.structure class to facilitate calculation of distances
    xs = xray.structure(crystal_symmetry=self.crystal_symmetry)
    for i, site in enumerate(self.sites):
      xs.add_scatterer(xray.scatterer("C%i" %i, site=site))

    xs = xs.sites_mod_short()
    sites_cart = xs.sites_cart()
    lengths = flex.double([matrix.col(sc).length() for sc in sites_cart])
    xs = xs.select(flex.sort_permutation(lengths))
    with open('peaks.pdb', 'wb') as f:
      print >> f, xs.as_pdb_file()

    vector_heights = flex.double()

    sites_frac = xs.sites_frac()
    pair_asu_table = xs.pair_asu_table(distance_cutoff=self.params.max_cell)
    asu_mappings = pair_asu_table.asu_mappings()
    distances = crystal.calculate_distances(pair_asu_table, sites_frac)
    vectors = []
    for di in distances:
      if di.distance < self.params.min_cell: continue
      i_seq, j_seq = di.i_seq, di.j_seq
      if i_seq > 0:
        continue
      # Is this the peak centred at (0,0,0)?
      assert (sum(x**2 for x in sites_frac[i_seq]) < 1e-8)
      rt_mx_ji = di.rt_mx_ji
      site_frac_ji = rt_mx_ji * sites_frac[j_seq]
      site_cart = self.fft_cell.orthogonalize(site_frac_ji)
      vectors.append(matrix.col(site_cart))
      #vector_heights.append(heights[j_seq])

    # XXX loop over these vectors and sort into groups similar to further down
    # group similar angle and lengths, also catch integer multiples of vectors

    vector_groups = []
    relative_length_tolerance = 0.1
    angle_tolerance = 10 # degrees

    orth = self.fft_cell.orthogonalize
    for v in vectors:
      length = v.length()
      if length < self.params.min_cell or length > self.params.max_cell:
        continue
      matched_group = False
      for group in vector_groups:
        mean_v = group.mean()
        mean_v_length = mean_v.length()
        if (abs(mean_v_length - length)/max(mean_v_length, length)
            < relative_length_tolerance):
          angle = mean_v.angle(v, deg=True)
          if angle < angle_tolerance:
            group.append(v, length)
            matched_group = True
            break
          elif abs(180-angle) < angle_tolerance:
            group.append(-v, length)
            matched_group = True
            break
      if not matched_group:
        group = vector_group()
        group.append(v, length)
        vector_groups.append(group)

    vectors = [g.mean() for g in vector_groups]

    # sort by length
    lengths = flex.double([v.length() for v in vectors])
    perm = flex.sort_permutation(lengths)
    vectors = [vectors[i] for i in perm]
    #vector_heights = [vector_heights[i] for i in perm]

    # exclude vectors that are (approximately) integer multiples of a shorter
    # vector
    unique_vectors = []
    for v in vectors:
      is_unique = True
      if i > 0:
        for v_u in unique_vectors:
          if is_approximate_integer_multiple(v_u, v):
            is_unique = False
            break
      if is_unique:
        unique_vectors.append(v)
    vectors = unique_vectors

    # choose the shortest vectors
    vectors = vectors[:30]

    self.candidate_basis_vectors = vectors

  def find_basis_vector_combinations_cluster_analysis(self):
    # hijack the xray.structure class to facilitate calculation of distances
    xs = xray.structure(crystal_symmetry=self.crystal_symmetry)
    for i, site in enumerate(self.sites):
      xs.add_scatterer(xray.scatterer("C%i" %i, site=site))

    xs = xs.sites_mod_short()
    xs = xs.select(xs.sites_frac().norms() < 0.45)
    cell_multiplier = 10
    xs1 = xs.customized_copy(
      unit_cell=uctbx.unit_cell([xs.unit_cell().parameters()[0]*cell_multiplier]*3))
    xs1.set_sites_cart(xs.sites_cart())
    xs = xs1
    sites_cart = xs.sites_cart()
    lengths = flex.double([matrix.col(sc).length() for sc in sites_cart])
    xs = xs.select(flex.sort_permutation(lengths))
    with open('peaks.pdb', 'wb') as f:
      print >> f, xs.as_pdb_file()

    vector_heights = flex.double()

    sites_frac = xs.sites_frac()
    pair_asu_table = xs.pair_asu_table(distance_cutoff=self.params.max_cell)
    asu_mappings = pair_asu_table.asu_mappings()
    distances = crystal.calculate_distances(pair_asu_table, sites_frac)
    vectors = []
    difference_vectors = []
    pairs = []
    for di in distances:
      if di.distance < self.params.min_cell: continue
      i_seq, j_seq = di.i_seq, di.j_seq
      if i_seq > j_seq: continue
      pairs.append((i_seq, j_seq))
      rt_mx_ji = di.rt_mx_ji
      site_frac_ji = rt_mx_ji * sites_frac[j_seq]
      site_cart_ji = xs.unit_cell().orthogonalize(site_frac_ji)
      site_cart_i = xs.unit_cell().orthogonalize(sites_frac[i_seq])
      vectors.append(matrix.col(site_cart_ji))
      diff_vec = matrix.col(site_cart_i) - matrix.col(site_cart_ji)
      if diff_vec[0] < 0:
        # only one hemisphere of difference vector space
        diff_vec = -diff_vec
      difference_vectors.append(diff_vec)

    if self.params.cluster_analysis.method == 'dbscan':
      i_cluster = self.cluster_analysis_dbscan(difference_vectors)
      min_cluster_size = 1
    elif self.params.cluster_analysis.method == 'hcluster':
      i_cluster = self.cluster_analysis_hcluster(difference_vectors)
      i_cluster -= 1 # hcluster starts counting at 1
      min_cluster_size = self.params.cluster_analysis.min_cluster_size

    if self.params.debug_plots:
      self.debug_plot_clusters(
        difference_vectors, i_cluster, min_cluster_size=min_cluster_size)


    clusters = []
    min_cluster_size = self.params.cluster_analysis.min_cluster_size
    for i in range(max(i_cluster)+1):
      isel = (i_cluster == i).iselection()
      if len(isel) < min_cluster_size:
        continue
      clusters.append(isel)

    cluster_point_sets = []
    centroids = []
    cluster_sizes = flex.int()

    difference_vectors = flex.vec3_double(difference_vectors)

    from libtbx.utils import flat_list
    for cluster in clusters:
      points = flat_list([pairs[i] for i in cluster])
      cluster_point_sets.append(set(points))
      d_vectors = difference_vectors.select(cluster)
      cluster_sizes.append(len(d_vectors))
      centroids.append(d_vectors.mean())

    # build a graph where each node is a centroid from the difference vector
    # cluster analysis above, and an edge is defined when there is a
    # significant overlap between the sets of peaks in the FFT map that
    # contributed to the difference vectors in two clusters
    import networkx as nx
    G = nx.Graph()
    G.add_nodes_from(range(len(cluster_point_sets)))

    cutoff_frac = 0.25
    for i in range(len(cluster_point_sets)):
      for j in range(i+1, len(cluster_point_sets)):
        intersection_ij = cluster_point_sets[i].intersection(
            cluster_point_sets[j])
        union_ij = cluster_point_sets[i].union(cluster_point_sets[j])
        frac_connected = len(intersection_ij)/len(union_ij)
        if frac_connected > cutoff_frac:
          G.add_edge(i, j)

    # iteratively find the maximum cliques in the graph
    # break from the loop if there are no cliques remaining or there are
    # fewer than 3 vectors in the remaining maximum clique
    # Allow 1 basis vector to be shared between two cliques, to allow for
    # cases where two lattices share one basis vectors (e.g. two plate
    # crystals exactly aligned in one direction, but not in the other two)
    distinct_cliques = []
    cliques = list(nx.find_cliques(G))
    cliques = sorted(cliques, key=len, reverse=True)
    for i, clique in enumerate(cliques):
      clique = set(clique)
      if len(clique) < 3:
        break
      is_distinct = True
      for c in distinct_cliques:
        if len(c.intersection(clique)) > 1:
          is_distinct = False
          break
      if is_distinct:
        distinct_cliques.append(clique)
        this_set = set()
        for i_cluster in clique:
          this_set = this_set.union(cluster_point_sets[i_cluster])
        print "Clique %i: %i lattice points" %(i+1, len(this_set))

    assert len(distinct_cliques) > 0

    print "Estimated number of lattices: %i" %len(distinct_cliques)

    self.candidate_basis_vectors = []
    self.candidate_crystal_models = []

    for clique in distinct_cliques:
      sel = flex.size_t(list(clique))
      vectors = flex.vec3_double(centroids).select(sel)
      perm = flex.sort_permutation(vectors.norms())
      vectors = [matrix.col(vectors[p]) for p in perm]

      # exclude vectors that are (approximately) integer multiples of a shorter
      # vector
      unique_vectors = []
      for v in vectors:
        is_unique = True
        for v_u in unique_vectors:
          if is_approximate_integer_multiple(v_u, v,
                                             relative_tolerance=0.01,
                                             angular_tolerance=0.5):
            is_unique = False
            break
        if is_unique:
          unique_vectors.append(v)
      vectors = unique_vectors

      self.candidate_basis_vectors.extend(vectors)
      candidate_orientation_matrices \
        = self.find_candidate_orientation_matrices(vectors, return_first=True)
      if len(candidate_orientation_matrices) == 0:
        continue
      # only take the first one
      crystal_model = candidate_orientation_matrices[0]
      # map to minimum reduced cell
      crystal_symmetry = crystal.symmetry(
        unit_cell=crystal_model.get_unit_cell(),
        space_group=crystal_model.get_space_group())
      cb_op = crystal_symmetry.change_of_basis_op_to_minimum_cell()
      crystal_model = crystal_model.change_basis(cb_op)
      self.candidate_crystal_models.append(crystal_model)

    for i_lattice in range(len(self.candidate_crystal_models)):
      if self.target_symmetry_primitive is not None:
        symmetrized_model = self.apply_symmetry(
          self.candidate_crystal_models[i_lattice], self.target_symmetry_primitive)
        print symmetrized_model.get_unit_cell()
        self.candidate_crystal_models[i_lattice] = symmetrized_model

    if self.params.debug:
      file_name = "vectors.pdb"
      a = self.params.max_cell
      cs = crystal.symmetry(unit_cell=(a,a,a,90,90,90), space_group="P1")
      xs = xray.structure(crystal_symmetry=cs)
      for v in difference_vectors:
        v = matrix.col(v)
        xs.add_scatterer(xray.scatterer("C", site=v/(a/10)))
      xs.sites_mod_short()
      with open(file_name, 'wb') as f:
        print >> f, xs.as_pdb_file()

    if self.params.debug:
      for crystal_model in self.candidate_crystal_models:
        print crystal_model

  def cluster_analysis_hcluster(self, vectors):
    from hcluster import linkage, fcluster
    import numpy

    params = self.params.cluster_analysis.hcluster
    X = numpy.array(vectors)
    linkage_method = params.linkage.method
    linkage_metric = params.linkage.metric
    criterion = params.cutoff_criterion
    Z = linkage(X, method=linkage_method, metric=linkage_metric)
    cutoff = params.cutoff
    i_cluster = fcluster(Z, cutoff, criterion=criterion)
    i_cluster = flex.int(i_cluster.astype(numpy.int32))
    return i_cluster

  def cluster_analysis_dbscan(self, vectors):
    import numpy as np

    from sklearn.cluster import DBSCAN
    from sklearn.preprocessing import StandardScaler

    vectors = flex.vec3_double(vectors)

    X = np.array(vectors)
    # scale the data - is this necessary/does it help or hinder?
    X = StandardScaler().fit_transform(X)

    # Compute DBSCAN
    params = self.params.cluster_analysis.dbscan
    db = DBSCAN(eps=params.eps, min_samples=params.min_samples).fit(X)
    core_samples = db.core_sample_indices_
    # core_samples is a list of numpy.int64 objects
    core_samples = flex.int([int(i) for i in core_samples])
    labels = flex.int(db.labels_.astype(np.int32))

    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

    print('Estimated number of clusters: %d' % n_clusters_)

    return labels

  def debug_plot_clusters(self, vectors, labels, min_cluster_size=1):
    assert len(vectors) == len(labels)
    from matplotlib import pyplot
    from mpl_toolkits.mplot3d import Axes3D
    import numpy

    # Black removed and is used for noise instead.
    fig = pyplot.figure()
    ax = fig.add_subplot(111, projection='3d')
    # label clusters smaller than min_cluster_size as noise
    unique_labels = set(labels)
    for k in unique_labels:
      sel_k = labels == k
      if (sel_k).count(True) < min_cluster_size:
        labels = labels.set_selected(sel_k, -1)
    unique_labels = set(labels)
    n_clusters = len(unique_labels) - (1 if -1 in unique_labels else 0)
    colours = pyplot.cm.Spectral(numpy.linspace(0, 1, len(unique_labels)))
    #vectors = flex.vec3_double(vectors)
    ax.scatter([0],[0],[0], c='k', marker='+', s=50)
    for k, col in zip(unique_labels, colours):
      class_members = (labels == k).iselection()
      if k == -1:# or len(class_members) < min_cluster_size:
          # Black used for noise.
        col = 'k'
        col = '0.25' # mid-grey
        markersize = 1
        marker = '+'
        #continue
      else:
        markersize = 10
        marker = 'o'
      if not isinstance(col, basestring) and len(col) == 4:
        # darken the edges
        frac = 0.75
        edgecolor = [col[0]*frac, col[1]*frac, col[2]*frac, col[3]]
      else:
        edgecolor = col
      vs = numpy.array([vectors[i] for i in class_members])
      xs = vs[:,0]
      ys = vs[:,1]
      zs = vs[:,2]
      # plot whole sphere for visual effect
      xs = numpy.concatenate((xs, -xs))
      ys = numpy.concatenate((ys, -ys))
      zs = numpy.concatenate((zs, -zs))
      ax.scatter(xs, ys, zs, c=col, marker=marker, s=markersize, edgecolor=edgecolor)

    pyplot.title('Estimated number of clusters: %d' % n_clusters)
    pyplot.show()
