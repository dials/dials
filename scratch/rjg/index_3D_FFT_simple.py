from __future__ import division
import sys
import cPickle as pickle
import math

try:
  # try importing scipy.linalg before any cctbx modules, otherwise we
  # sometimes get a segmentation fault/core dump if it is imported after
  # scipy.linalg is a dependency of sklearn.cluster.DBSCAN
  import scipy.linalg # import dependency
except ImportError, e:
  pass

import iotbx.phil
from scitbx import matrix
from scitbx import fftpack

from cctbx.array_family import flex
from cctbx import crystal, miller, sgtbx, uctbx, xray

import dxtbx
from dials.algorithms.centroid import centroid_px_to_mm
from dials.model.data import ReflectionList
from dials.model.experiment.crystal_model import Crystal

import libtbx.load_env
dials_path = libtbx.env.dist_path('dials')

master_phil_scope = iotbx.phil.parse("""
min_cell = 20
  .type = float(value_min=0)
  .help = "Minimum length of candidate unit cell basis vectors (in Angstrom)."
max_cell = 160
  .type = float(value_min=0)
  .help = "Maximum length of candidate unit cell basis vectors (in Angstrom)."
reciprocal_space_grid {
  n_points = 256
    .type = int(value_min=0)
  d_min = 4
    .type = float(value_min=0)
    .help = "The high resolution limit in Angstrom for spots to include in "
            "the initial indexing."
}
b_iso = 200
  .type = float(value_min=0)
rmsd_cutoff = 15
  .type = float(value_min=0)
scan_range = None
  .help = "The range of images to use in indexing. Number of arguments"
          "must be a factor of two. Specifying \"0 0\" will use all images"
          "by default. The given range follows C conventions"
          "(e.g. j0 <= j < j1)."
  .type = ints(size=2)
  .multiple = True
known_symmetry {
  space_group = None
    .type = space_group
  unit_cell = None
    .type = unit_cell
  relative_length_tolerance = 0.1
    .type = float
    .help = "Relative tolerance for unit cell lengths in unit cell comparision."
  absolute_angle_tolerance = 10
    .type = float
    .help = "Angular tolerance (in degrees) in unit cell comparison."
}
debug = False
  .type = bool
debug_plots = False
  .type = bool
  .help = "Requires matplotlib"
include file %s/data/refinement.phil
refinement_protocol {
  n_macro_cycles = 3
    .type = int(value_min=0)
  d_min_step = 0.5
    .type = float(value_min=0.0)
    .help = "Reduction per step in d_min for reflections to include in refinement."
  d_min_final = None
    .type = float(value_min=0.0)
    .help = "Do not ever include reflections below this value in refinement."
  verbosity = 1
    .type = int(value_min=0)
}
export_xds_files = False
  .type = bool
  .help = "Export results as XDS.INP, XPARM.XDS for integration with XDS."
multiple_lattice_search = False
  .type = bool
cluster_analysis {
  method = *dbscan hcluster
    .type = choice
  hcluster {
    linkage {
      method = *ward
        .type = choice
      metric = *euclidean
        .type = choice
    }
    cutoff = 15
      .type = float(value_min=0)
    cutoff_criterion = *distance inconsistent
      .type = choice
  }
  dbscan {
    eps = 0.05
      .type = float(value_min=0.0)
    min_samples = 30
      .type = int(value_min=1)
  }
  min_cluster_size = 20
    .type = int(value_min=0)
  intersection_union_ratio_cutoff = 0.4
    .type = float(value_min=0.0, value_max=1.0)
}
""" %dials_path, process_includes=True)

master_params = master_phil_scope.fetch().extract()


class vector_group(object):
  def __init__(self):
    self.vectors = []
    self.lengths = []
    self._mean = None

  def append(self, vector, length):
    self.vectors.append(vector)
    self.lengths.append(length)
    self._mean = self.compute_mean()

  def mean(self):
    if self._mean is None:
      self._mean = self.compute_mean()
    return self._mean

  def compute_mean(self):
    sum_x = 0
    sum_y = 0
    sum_z = 0
    for v in self.vectors:
      sum_x += v.elems[0]
      sum_y += v.elems[1]
      sum_z += v.elems[2]
    return matrix.col((sum_x, sum_y, sum_z))/len(self.vectors)


def is_approximate_integer_multiple(vec_a, vec_b,
                                     relative_tolerance=0.2,
                                     angular_tolerance=5.0):
  length_a = vec_a.length()
  length_b = vec_b.length()
  assert length_b >= length_a
  angle = vec_a.angle(vec_b, deg=True)
  if angle < angular_tolerance or abs(180-angle) < angular_tolerance:
    n = length_b/length_a
    if abs(round(n) - n) < relative_tolerance:
      return True
  return False


deg_to_radians = math.pi/180


class indexer(object):

  def __init__(self, reflections, sweep, params=None):
    self.reflections = reflections
    # the lattice a given reflection belongs to: a value of -1 indicates
    # that a reflection doesn't belong to any lattice so far
    self.reflections_i_lattice = flex.int(self.reflections.size(), -1)
    self.sweep = sweep
    self.goniometer = sweep.get_goniometer()
    self.detector = sweep.get_detector()
    self.scan = sweep.get_scan()
    self.beam = sweep.get_beam()
    if params is None: params = master_params
    self.params = params

    self.target_symmetry_primitive = None
    self.target_symmetry_centred = None
    if (self.params.known_symmetry.space_group is not None or
        self.params.known_symmetry.unit_cell is not None):
      is_centred = False
      if self.params.known_symmetry.space_group is not None:
        space_group_info = self.params.known_symmetry.space_group
        is_centred = space_group_info.group().conventional_centring_type_symbol() != 'P'
        cb_op_to_primitive = space_group_info.change_of_basis_op_to_primitive_setting()
        sgi_primitive = space_group_info.change_basis(cb_op_to_primitive)
        space_group_primitive = sgi_primitive.group()
      else:
        space_group_primitive = sgtbx.space_group("P 1")
      self.target_symmetry_primitive = crystal.symmetry(
        unit_cell=self.params.known_symmetry.unit_cell,
        space_group=space_group_primitive,
        assert_is_compatible_unit_cell=False)
      if is_centred and self.params.known_symmetry.unit_cell is not None:
        self.target_symmetry_centred = self.target_symmetry_primitive.change_basis(
          self.target_symmetry_primitive.change_of_basis_op_to_reference_setting())
      if self.params.known_symmetry.unit_cell is not None:
        assert self.target_symmetry_primitive.unit_cell().is_similar_to(
          self.params.known_symmetry.unit_cell)

  def index(self):
    n_points = self.params.reciprocal_space_grid.n_points
    self.gridding = fftpack.adjust_gridding_triple(
      (n_points,n_points,n_points), max_prime=5)
    n_points = self.gridding[0]
    self.prepare_reflections()
    self.filter_reflections_by_scan_range()
    self.map_centroids_to_reciprocal_space_grid()

    print "Number of centroids used: %i" %(
      (self.reciprocal_space_grid>0).count(True))
    self.fft()
    if self.params.debug:
      self.debug_write_reciprocal_lattice_points_as_pdb()
      self.debug_write_ccp4_map(map_data=self.grid_real, file_name="patt.map")
    self.find_peaks()
    if self.params.multiple_lattice_search:
      self.find_basis_vector_combinations_cluster_analysis()
      if self.params.debug:
        self.debug_show_candidate_basis_vectors()
      crystal_models = self.candidate_crystal_models
    else:
      self.find_candidate_basis_vectors()
      if self.params.debug:
        self.debug_show_candidate_basis_vectors()
      self.candidate_crystal_models = self.find_candidate_orientation_matrices(
        self.candidate_basis_vectors)
      crystal_models = self.candidate_crystal_models[:1]

    self.refined_crystal_models = []
    for i_lattice, crystal_model in enumerate(crystal_models):
      self.i_lattice = i_lattice

      print
      print "#" * 80
      print "Starting refinement of crystal model %i" %(i_lattice+1)
      print "Starting crystal model:"
      print crystal_model
      print "#" * 80

      self.d_min = self.params.reciprocal_space_grid.d_min
      self.indexed_reflections, _ = self.index_reflections_given_orientation_matix(
        crystal_model, verbose=1)
      if self.target_symmetry_primitive is not None:
        symmetrized_model = self.apply_symmetry(
          crystal_model, self.target_symmetry_primitive)
        self.indexed_reflections, _ = self.index_reflections_given_orientation_matix(
          symmetrized_model, verbose=1)
        crystal_model = symmetrized_model
      for i_cycle in range(self.params.refinement_protocol.n_macro_cycles):
        print "Starting refinement (macro-cycle %i)" %(i_cycle+1)
        print
        self.refine(crystal_model)

        self.d_min -= self.params.refinement_protocol.d_min_step
        self.d_min = max(self.d_min, self.params.refinement_protocol.d_min_final)
        print "Increasing resolution to %.1f Angstrom" %self.d_min
        n_indexed_last_cycle = self.indexed_reflections.size()

        self.indexed_reflections, _ = self.index_reflections_given_orientation_matix(
          crystal_model, verbose=1)
        # reset in case some reflections are no longer indexed by this lattice
        self.reflections_i_lattice.set_selected(
          self.reflections_i_lattice == i_lattice, -1)
        self.reflections_i_lattice.set_selected(
          self.indexed_reflections, i_lattice)
        if self.indexed_reflections.size() == n_indexed_last_cycle:
          print "No more reflections indexed this cycle - finished with refinement"
          break
        elif self.d_min == self.params.refinement_protocol.d_min_final:
          print "Target d_min_final reached: finished with refinement"
          break

      print "Final crystal model:"
      print crystal_model

      self.refined_crystal_models.append(crystal_model)
      suffix = ""
      if len(crystal_models) > 1:
        suffix = "_%i" %(i_lattice+1)
      self.export_as_json(crystal_model, self.sweep, suffix=suffix)
      if self.params.export_xds_files:
        self.export_xds_files(crystal_model, self.sweep, suffix=suffix)
      self.export_reflections(file_name='indexed%s.pickle' %suffix,
                              indexed_only=True)

    #self.export_reflections(indexed_only=False)
    if self.params.debug:
      self.predict_reflections(self.candidate_crystal_models[0])
      self.export_predicted_reflections()

    print "Final refined crystal models:"
    for i, crystal_model in enumerate(self.refined_crystal_models):
      print "model %i (%i reflections):" %(
        i+1, (self.reflections_i_lattice == i).count(True))
      print crystal_model

  def prepare_reflections(self):
    """Reflections that come from dials.spotfinder only have the centroid
    position and variance set, """

    for i_ref, refl in enumerate(self.reflections):
      # just a quick check for now that the reflections haven't come from
      # somewhere else
      assert refl.image_coord_mm == (0,0)

      # set reflection properties that might be needed by the dials refinement
      # engine, and convert values from pixels and image number to mm/rads
      refl.frame_number = refl.centroid_position[2]
      centroid_position, centroid_variance, _ = centroid_px_to_mm(
        self.detector, self.scan,
        refl.centroid_position,
        refl.centroid_variance,
        (1,1,1))
      refl.centroid_position = centroid_position
      refl.centroid_variance = centroid_variance
      refl.rotation_angle = centroid_position[2]

  def filter_reflections_by_scan_range(self):
    self.reflections_in_scan_range = flex.size_t()
    for i_ref, refl in enumerate(self.reflections):
      frame_number = refl.frame_number

      if len(self.params.scan_range):
        reflections_in_range = False
        for scan_range in self.params.scan_range:
          if scan_range is None: continue
          range_start, range_end = scan_range
          if frame_number >= range_start and frame_number < range_end:
            reflections_in_range = True
            break
        if not reflections_in_range:
          continue
      self.reflections_in_scan_range.append(i_ref)

  def map_centroids_to_reciprocal_space(self):
    self.reciprocal_space_points = flex.vec3_double()

    wavelength = self.beam.get_wavelength()
    s0 = matrix.col(self.beam.get_s0())
    rotation_axis = matrix.col(
      self.goniometer.get_rotation_axis())

    for i_ref in self.reflections_in_scan_range:
      refl = self.reflections[i_ref]
      s1 = matrix.col(
        self.detector.get_lab_coord(refl.centroid_position[:2]))
      s1 = s1.normalize()/wavelength
      refl.beam_vector = tuple(s1) # needed by ray_intersection
      S = s1 - s0
      phi = refl.rotation_angle
      point = S.rotate_around_origin(rotation_axis, -phi, deg=False)
      self.reciprocal_space_points.append(tuple(point))

  def map_centroids_to_reciprocal_space_grid(self):
    self.map_centroids_to_reciprocal_space()
    assert len(self.reciprocal_space_points) == len(self.reflections_in_scan_range)
    wavelength = self.beam.get_wavelength()
    d_min = self.params.reciprocal_space_grid.d_min

    n_points = self.gridding[0]
    rlgrid = 2 / (d_min * n_points)

    # real space FFT grid dimensions
    cell_lengths = [n_points * d_min/2 for i in range(3)]
    self.unit_cell = uctbx.unit_cell(cell_lengths+[90]*3)
    self.crystal_symmetry = crystal.symmetry(unit_cell=self.unit_cell,
                                             space_group_symbol="P1")

    print "FFT gridding: (%i,%i,%i)" %self.gridding

    grid = flex.double(flex.grid(self.gridding), 0)

    reflections_used_for_indexing = flex.size_t()

    for i_pnt, point in enumerate(self.reciprocal_space_points):
      point = matrix.col(point)
      spot_resolution = 1/point.length()
      if spot_resolution < d_min:
        continue

      grid_coordinates = [int(round(point[i]/rlgrid)+n_points/2) for i in range(3)]
      if max(grid_coordinates) >= n_points: continue # this reflection is outside the grid
      if min(grid_coordinates) < 0: continue # this reflection is outside the grid
      T = math.exp(-self.params.b_iso * point.length()**2 / 4)
      grid[grid_coordinates] = T

      i_ref = self.reflections_in_scan_range[i_pnt]
      reflections_used_for_indexing.append(i_ref)

    self.reciprocal_space_grid = grid
    self.reflections_used_for_indexing = reflections_used_for_indexing

  def fft(self):
    fft = fftpack.complex_to_complex_3d(self.gridding)
    grid_complex = flex.complex_double(
      reals=self.reciprocal_space_grid,
      imags=flex.double(self.reciprocal_space_grid.size(), 0))
    self.grid_transformed = fft.forward(grid_complex)
    #self.grid_real = flex.pow2(flex.abs(self.grid_transformed))
    self.grid_real = flex.pow2(flex.real(self.grid_transformed))
    #self.grid_real = flex.pow2(flex.imag(self.grid_transformed))

  def find_peaks(self):
    grid_real = self.grid_real.deep_copy()
    rmsd = math.sqrt(
      flex.mean(flex.pow2(grid_real.as_1d()-flex.mean(grid_real.as_1d()))))
    grid_real.set_selected(grid_real < (self.params.rmsd_cutoff)*rmsd, 0)
    grid_real_binary = grid_real.deep_copy()
    grid_real_binary.as_1d().set_selected(grid_real.as_1d() > 0, 1)
    grid_real_binary = grid_real_binary.iround()
    from cctbx import masks
    flood_fill = masks.flood_fill(grid_real_binary, self.unit_cell)
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

    self.sites = flood_fill.centres_of_mass_frac().select(isel)

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
      site_cart = self.unit_cell.orthogonalize(site_frac_ji)
      vectors.append(matrix.col(site_cart))
      #vector_heights.append(heights[j_seq])

    # XXX loop over these vectors and sort into groups similar to further down
    # group similar angle and lengths, also catch integer multiples of vectors

    vector_groups = []
    relative_length_tolerance = 0.1
    angle_tolerance = 10 # degrees

    orth = self.unit_cell.orthogonalize
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
      site_cart_ji = self.unit_cell.orthogonalize(site_frac_ji)
      site_cart_i = self.unit_cell.orthogonalize(sites_frac[i_seq])
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

    difference_vectors = flex.vec3_double(difference_vectors)

    from libtbx.utils import flat_list
    for cluster in clusters:
      points = flat_list([pairs[i] for i in cluster])
      cluster_point_sets.append(set(points))
      centroids.append(difference_vectors.select(cluster).mean())

    combinations = flex.vec3_int()
    unions = []
    intersections = []
    n_common_points = flex.int()
    fraction_common_points = flex.double()

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
    for clique in cliques:
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

    assert len(distinct_cliques) > 0

    print "Estimated number of lattices: %i" %len(distinct_cliques)

    self.candidate_basis_vectors = []
    self.candidate_crystal_models = []

    for clique in distinct_cliques:
      vectors = flex.vec3_double(centroids).select(flex.size_t(list(clique)))
      perm = flex.sort_permutation(vectors.norms())
      vectors = [matrix.col(vectors[p]) for p in perm]
      self.candidate_basis_vectors.extend(vectors)
      candidate_orientation_matrices \
        = self.find_candidate_orientation_matrices(vectors, return_first=True)
      # only take the first one
      crystal_model = candidate_orientation_matrices[0]
      # map to minimum reduced cell
      crystal_symmetry = crystal.symmetry(
        unit_cell=crystal_model.get_unit_cell(),
        space_group=crystal_model.get_space_group())
      cb_op = crystal_symmetry.change_of_basis_op_to_minimum_cell()
      crystal_model = crystal_model.change_basis(cb_op)
      self.candidate_crystal_models.append(crystal_model)

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

  def cluster_analysis_hcluster(self, vectors):
    from hcluster import pdist, linkage, dendrogram, fcluster
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

  def find_candidate_orientation_matrices(self, candidate_basis_vectors,
                                          return_first=False):
    candidate_crystal_models = []
    vectors = candidate_basis_vectors

    min_angle = 20 # degrees, aritrary cutoff
    for i in range(len(vectors)):
      a = vectors[i]
      for j in range(i, len(vectors)):
        b = vectors[j]
        angle = a.angle(b, deg=True)
        if angle < min_angle or (180-angle) < min_angle:
          continue
        a_cross_b = a.cross(b)
        gamma = a.angle(b, deg=True)
        if gamma < 90:
          # all angles obtuse if possible please
          b = -b
          gamma = 180 - gamma
          a_cross_b = -a_cross_b
        for k in range(j, len(vectors)):
          c = vectors[k]
          if abs(90-a_cross_b.angle(c, deg=True)) < min_angle:
            continue
          alpha = b.angle(c, deg=True)
          if alpha < 90:
            c = -c
          beta = c.angle(a, deg=True)
          if a_cross_b.dot(c) < 0:
            # we want right-handed basis set, therefore invert all vectors
            a = -a
            b = -b
            c = -c
            #assert a.cross(b).dot(c) > 0
          model = Crystal(a, b, c, space_group_symbol="P 1")
          uc = model.get_unit_cell()
          if self.target_symmetry_primitive is not None:
            symmetrized_model = self.apply_symmetry(
              model, self.target_symmetry_primitive)
            if symmetrized_model is None:
              if self.target_symmetry_centred is not None:
                symmetrized_model = self.apply_symmetry(
                model, self.target_symmetry_centred,
                return_primitive_setting=True)
                if symmetrized_model is not None:
                  model = symmetrized_model
                  uc = model.get_unit_cell()
              if symmetrized_model is None:
                continue

          params = uc.parameters()
          if uc.volume() > (params[0]*params[1]*params[2]/100):
            # unit cell volume cutoff from labelit 2004 paper
            candidate_crystal_models.append(model)
            if return_first:
              return candidate_crystal_models
    return candidate_crystal_models

  def predict_reflections(self, crystal_model):
    from dials.algorithms.integration import ReflectionPredictor
    predictor = ReflectionPredictor()

    sigma_divergence = self.beam.get_sigma_divergence()
    mosaicity = crystal_model.get_mosaicity()

    if sigma_divergence == 0.0:
      self.beam.set_sigma_divergence(0.02) # degrees
    if mosaicity == 0.0:
      crystal_model.set_mosaicity(0.139) # degrees

    reflections = predictor(self.sweep, crystal_model)
    self.predicted_reflections = reflections
    return self.predicted_reflections

  def export_predicted_reflections(self, file_name='predictions.pickle'):
    from dials.model.serialize import dump
    dump.reflections(self.predicted_reflections, file_name)

  def apply_symmetry(self, crystal_model, target_symmetry,
                     return_primitive_setting=False):
    unit_cell = crystal_model.get_unit_cell()
    target_unit_cell = target_symmetry.unit_cell()
    target_space_group = target_symmetry.space_group()
    A = crystal_model.get_A()
    A_inv = A.inverse()
    unit_cell_is_similar = False
    real_space_a = A_inv[:3]
    real_space_b = A_inv[3:6]
    real_space_c = A_inv[6:9]
    basis_vectors = [real_space_a, real_space_b, real_space_c]
    min_bmsd = 1e8
    best_perm = None
    for perm in ((0,1,2), (1,2,0), (2,0,1)):
      crystal_model = Crystal(
        basis_vectors[perm[0]],
        basis_vectors[perm[1]],
        basis_vectors[perm[2]],
        space_group=target_space_group)
      unit_cell = crystal_model.get_unit_cell()
      uc = target_unit_cell
      if uc is None:
        uc = unit_cell
      # XXX what about permuting the target_unit_cell (if not None)?
      symm_target_sg = crystal.symmetry(
        unit_cell=uc,
        space_group=target_space_group,
        assert_is_compatible_unit_cell=False)
      # this assumes that the initial basis vectors are good enough that
      # we can tell which should be the unique axis - probably not a robust
      # solution
      if unit_cell.is_similar_to(
        symm_target_sg.unit_cell(),
        relative_length_tolerance=self.params.known_symmetry.relative_length_tolerance,
        absolute_angle_tolerance=self.params.known_symmetry.absolute_angle_tolerance):
        bmsd = unit_cell.bases_mean_square_difference(
          symm_target_sg.unit_cell())
        min_bmsd = min(min_bmsd, bmsd)
        if min_bmsd == bmsd:
          best_perm = list(perm)
    if best_perm is None:
      return None
    crystal_model = Crystal(
      basis_vectors[best_perm[0]],
      basis_vectors[best_perm[1]],
      basis_vectors[best_perm[2]],
      space_group=target_space_group)
    model = crystal_model
    #cb_op_target_ref = symm_target_sg.space_group_info().type().cb_op()
    #symm_target_sg_ref = symm_target_sg.change_basis(cb_op_target_ref)
    from rstbx.symmetry.constraints import parameter_reduction
    s = parameter_reduction.symmetrize_reduce_enlarge(target_space_group)
    s.set_orientation(crystal_model.get_A())
    #s.symmetrize()
    #direct_matrix = s.orientation.change_basis(cb_op_target_ref).direct_matrix()
    if return_primitive_setting:
      sgi = sgtbx.space_group_info(group=target_space_group)
      cb_op_to_primitive = sgi.change_of_basis_op_to_primitive_setting()
      direct_matrix = s.orientation.change_basis(
        cb_op_to_primitive).direct_matrix()
    else:
      direct_matrix = s.orientation.direct_matrix()
    a = matrix.col(direct_matrix[:3])
    b = matrix.col(direct_matrix[3:6])
    c = matrix.col(direct_matrix[6:9])
    ## verify it is still right-handed basis set
    #assert a.cross(b).dot(c) > 0
    model = Crystal(
      a, b, c, space_group=target_space_group)
    return model

  def index_reflections_given_orientation_matix(
      self, crystal_model, tolerance=0.2, verbose=0, plot_differences=False):

    if verbose > 1:
      print "Candidate crystal model:"
      print crystal_model

    n_rejects = 0

    miller_indices = flex.miller_index()
    indexed_reflections = flex.size_t()

    A = crystal_model.get_A()
    A_inv = A.inverse()

    diff_h = flex.double()
    diff_k = flex.double()
    diff_l = flex.double()

    for i_ref in self.reflections_in_scan_range:
      if self.reflections_i_lattice[i_ref] > -1:
        # this reflection has already been indexed by a previous lattice
        continue
      i_rlp = flex.first_index(self.reflections_in_scan_range, i_ref)
      rlp = self.reciprocal_space_points[i_rlp]
      rlp = matrix.col(rlp)
      spot_resolution = 1/rlp.length()
      if spot_resolution < self.d_min:
        continue
      refl = self.reflections[i_ref]
      hkl_float = A_inv * rlp
      hkl_int = [int(round(h)) for h in hkl_float]
      if plot_differences:
        diff = matrix.col(hkl_int) - hkl_float
        diff_h.append(diff[0])
        diff_k.append(diff[1])
        diff_l.append(diff[2])
      max_difference = max([abs(hkl_float[i] - hkl_int[i]) for i in range(3)])
      if max_difference > tolerance:
        n_rejects += 1
        continue
      miller_indices.append(hkl_int)
      refl.miller_index = hkl_int
      indexed_reflections.append(i_ref)

    #n_rejects = n_rejects
    if verbose > 0:
      print "%i reflections indexed successfully (%i rejects)" %(
      indexed_reflections.size(), n_rejects)

    if plot_differences:
      from matplotlib import pyplot
      n_slots = 50
      hist_h = flex.histogram(diff_h, n_slots=n_slots)
      hist_k = flex.histogram(diff_k, n_slots=n_slots)
      hist_l = flex.histogram(diff_l, n_slots=n_slots)
      pyplot.plot(hist_h.slot_centers(), hist_h.slots())
      pyplot.plot(hist_k.slot_centers(), hist_k.slots())
      pyplot.plot(hist_l.slot_centers(), hist_l.slots())
      pyplot.show()

    return indexed_reflections, miller_indices

  def refine(self, crystal_model):
    from dials.algorithms.spot_prediction import ray_intersection
    reflections_for_refinement = ray_intersection(
      self.detector, self.reflections.select(self.indexed_reflections))
    verbosity = self.params.refinement_protocol.verbosity

    params = self.params.refinement
    from dials.algorithms.refinement import RefinerFactory
    refine = RefinerFactory.from_parameters(self.params, verbosity)
    scan_range_min = max(
      int(math.floor(flex.min(self.reflections.select(
        self.reflections_in_scan_range).frame_number()))),
      self.sweep.get_array_range()[0])
    scan_range_max = min(
      int(math.ceil(flex.max(self.reflections.select(
        self.reflections_in_scan_range).frame_number()))),
      self.sweep.get_array_range()[1])
    sweep = self.sweep[scan_range_min:scan_range_max]
    refine.prepare(sweep, crystal_model, reflections_for_refinement)
    #rmsds = refine.rmsds()
    refined = refine()

    if not (params.parameterisation.beam.fix_beam
            and params.parameterisation.detector.fix_detector):
      # Experimental geometry may have changed - re-map centroids to
      # reciprocal space
      self.map_centroids_to_reciprocal_space()

  def debug_show_candidate_basis_vectors(self):

    vectors = self.candidate_basis_vectors

    for i, v in enumerate(vectors):
      print i, v.length()# , vector_heights[i]

    # print a table of the angles between each pair of vectors

    angles = flex.double(len(vectors)**2)
    angles.reshape(flex.grid(len(vectors), len(vectors)))

    for i in range(len(vectors)):
      v_i = vectors[i]
      for j in range(i+1, len(vectors)):
        v_j = vectors[j]
        angles[i,j] = v_i.angle(v_j, deg=True)

    print (" "*7),
    for i in range(len(vectors)):
      print "%7.3f" % vectors[i].length(),
    print
    for i in range(len(vectors)):
      print "%7.3f" % vectors[i].length(),
      for j in range(len(vectors)):
        if j <= i:
          print (" "*7),
        else:
          print "%5.1f  " %angles[i,j],
      print

  def debug_write_reciprocal_lattice_points_as_pdb(
      self, file_name='reciprocal_lattice.pdb'):
    from cctbx import crystal, xray
    cs = crystal.symmetry(unit_cell=(1000,1000,1000,90,90,90), space_group="P1")
    xs = xray.structure(crystal_symmetry=cs)
    for site in self.reciprocal_space_points:
      xs.add_scatterer(xray.scatterer("C", site=site))

    xs.sites_mod_short()
    with open(file_name, 'wb') as f:
      print >> f, xs.as_pdb_file()

  def debug_write_ccp4_map(self, map_data, file_name):
    from iotbx import ccp4_map
    gridding_first = (0,0,0)
    gridding_last = map_data.all()
    labels = ["cctbx.miller.fft_map"]
    ccp4_map.write_ccp4_map(
      file_name=file_name,
      unit_cell=self.unit_cell,
      space_group=sgtbx.space_group("P1"),
      gridding_first=gridding_first,
      gridding_last=gridding_last,
      map_data=map_data,
      labels=flex.std_string(labels))

  def export_as_json(self, crystal_model, sweep, suffix=None, compact=False):
    from dials.model.serialize.dump import crystal as dump_crystal
    from dxtbx.serialize import dump
    if suffix is None:
      suffix = ''
    with open('crystal%s.json' %suffix, 'wb') as f:
      dump_crystal(crystal_model, f, compact=compact)
    with open('sweep%s.json' %suffix, 'wb') as f:
      dump.imageset(sweep, f, compact=compact)

  def export_reflections(self, file_name="indexed.pickle", indexed_only=False):
    reflections = self.reflections
    if indexed_only:
      reflections = reflections.select(self.indexed_reflections)
    with open(file_name, 'wb') as f:
      pickle.dump(reflections, f)

  def export_xds_files(self, crystal_model, sweep, suffix=None):
    from dxtbx.serialize import xds
    if suffix is None:
      suffix = ''
    crystal_model = crystal_model.change_basis(
      crystal_model.get_space_group().info().change_of_basis_op_to_reference_setting())
    A = crystal_model.get_A()
    A_inv = A.inverse()
    real_space_a = A_inv.elems[:3]
    real_space_b = A_inv.elems[3:6]
    real_space_c = A_inv.elems[6:9]
    to_xds = xds.to_xds(sweep)
    with open('XDS%s.INP' %suffix, 'wb') as f:
      to_xds.XDS_INP(out=f, job_card="XYCORR INIT DEFPIX INTEGRATE CORRECT")
    with open('XPARM%s.XDS' %suffix, 'wb') as f:
      to_xds.xparm_xds(
        real_space_a, real_space_b, real_space_c,
        crystal_model.get_space_group().type().number(),
        out=f)

def run(args):
  import time
  from libtbx.phil import command_line
  from dxtbx.imageset import ImageSetFactory

  args = sys.argv[1:]
  sweeps = ImageSetFactory.new(args, ignore_unknown=True)
  assert len(sweeps) == 1
  sweep = sweeps[0]
  sweep_filenames = sweep.paths()
  args = [arg for arg in args if arg not in sweep_filenames]

  cmd_line = command_line.argument_interpreter(master_params=master_phil_scope)
  working_phil, args = cmd_line.process_and_fetch(
      args=args, custom_processor="collect_remaining")
  working_phil.show()

  reflections_filename = args[0]

  gonio = sweep.get_goniometer()
  detector = sweep.get_detector()
  scan = sweep.get_scan()
  beam = sweep.get_beam()
  print detector
  print scan
  print gonio
  print beam

  t1 = time.time()
  with open(reflections_filename, 'rb') as f:
    reflections = pickle.load(f)
  t2 = time.time()
  print "Time taken loading reflection file: %.3fs" %(t2-t1)

  idxr = indexer(reflections, sweep,
                 params=working_phil.extract())
  idxr.index()
  return


if __name__ == '__main__':
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run(sys.argv[1:])
