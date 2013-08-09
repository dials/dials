from __future__ import division
import sys
import cPickle as pickle
import math

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
fft_n_points = 256
  .type = int(value_min=0)
d_min = 4
  .type = float(value_min=0)
  .help = "The high resolution limit in Angstrom for spots to include in indexing."
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
include file %s/data/refinement.phil
refinement {
  n_macro_cycles = 3
    .type = int(value_min=0)
  verbosity = 1
    .type = int(value_min=0)
}
export_xds_files = False
  .type = bool
  .help = "Export results as XDS.INP, XPARM.XDS for integration with XDS."
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
      if is_centred:
        self.target_symmetry_centred = self.target_symmetry_primitive.change_basis(
          self.target_symmetry_primitive.change_of_basis_op_to_reference_setting())
      if self.params.known_symmetry.unit_cell is not None:
        assert self.target_symmetry_primitive.unit_cell().is_similar_to(
          self.params.known_symmetry.unit_cell)

  def index(self):
    n_points = self.params.fft_n_points
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
    self.find_candidate_basis_vectors()
    if self.params.debug:
      self.debug_show_candidate_basis_vectors()
    self.find_candidate_orientation_matrices()
    crystal_model = self.candidate_crystal_models[0]
    self.d_min = self.params.d_min
    self.index_reflections_given_orientation_matix(
      crystal_model)
    if self.target_symmetry_primitive is not None:
      symmetrized_model = self.apply_symmetry(
        crystal_model, self.target_symmetry_primitive)
      self.index_reflections_given_orientation_matix(symmetrized_model)
      crystal_model = symmetrized_model
    for i in range(self.params.refinement.n_macro_cycles):
      print "Starting refinement (macro-cycle %i)" %(i+1)
      print
      self.refine(crystal_model)
      self.index_reflections_given_orientation_matix(crystal_model)

    self.candidate_crystal_models.insert(0, crystal_model)

    self.export_as_json()
    self.export_reflections(indexed_only=False)
    if self.params.debug:
      self.predict_reflections(self.candidate_crystal_models[0])
      self.export_predicted_reflections()
    if self.params.export_xds_files:
      self.export_xds_files()

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

    n_points = self.gridding[0]
    rlgrid = 2 / (self.params.d_min * n_points)

    # real space FFT grid dimensions
    cell_lengths = [n_points * self.params.d_min/2 for i in range(3)]
    self.unit_cell = uctbx.unit_cell(cell_lengths+[90]*3)
    self.crystal_symmetry = crystal.symmetry(unit_cell=self.unit_cell,
                                             space_group_symbol="P1")

    print "FFT gridding: (%i,%i,%i)" %self.gridding

    grid = flex.double(flex.grid(self.gridding), 0)

    reflections_used_for_indexing = flex.size_t()

    for i_pnt, point in enumerate(self.reciprocal_space_points):
      point = matrix.col(point)
      spot_resolution = 1/point.length()
      if spot_resolution < self.params.d_min:
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
    self.grid_real = flex.pow2(flex.abs(self.grid_transformed))

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
        0.3 * flex.max(flood_fill.grid_points_per_void()[1:]))).iselection()
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

  def find_candidate_orientation_matrices(self):
    self.candidate_crystal_models = []
    vectors = self.candidate_basis_vectors

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
              continue

          params = uc.parameters()
          if uc.volume() > (params[0]*params[1]*params[2]/100):
            # unit cell volume cutoff from labelit 2004 paper
            self.candidate_crystal_models.append(model)

  def predict_reflections(self, crystal_model):
    from dials.algorithms.spot_prediction import IndexGenerator
    from dials.algorithms.spot_prediction import ray_intersection
    from dials.algorithms.spot_prediction import reflection_frames
    from dials.algorithms.shoebox import BBoxCalculator
    from dials.algorithms.spot_prediction import RayPredictor
    from math import pi

    s0 = self.beam.get_s0()
    m2 = self.goniometer.get_rotation_axis()
    UB = crystal_model.get_U() * crystal_model.get_B()
    dphi = self.scan.get_oscillation_range(deg=False)

    d_min = self.detector.get_max_resolution_at_corners(
      s0, self.beam.get_wavelength())

    index_generator = IndexGenerator(
      crystal_model.get_unit_cell(),
      sgtbx.space_group("P 1").type(), d_min)
    miller_indices = index_generator.to_array()

    predict_rays = RayPredictor(s0, m2, dphi)
    self.predicted_reflections = reflection_frames(self.scan, ray_intersection(
        self.detector, predict_rays(miller_indices, UB)))

    sigma_divergence = self.beam.get_sigma_divergence()
    mosaicity = crystal_model.get_mosaicity()

    if sigma_divergence == 0.0:
      sigma_divergence = 0.02 # degrees
    if mosaicity == 0.0:
      mosaicity = 0.1 # degrees

    # Set the divergence and mosaicity
    n_sigma = 5.0
    delta_divergence = n_sigma * sigma_divergence * pi / 180.0
    delta_mosaicity = n_sigma * mosaicity * pi / 180.0

    # Create the bounding box calculator
    calculate_bbox = BBoxCalculator(self.beam, self.detector, self.goniometer,
        self.scan, delta_divergence, delta_mosaicity)

    # Calculate the frame numbers of all the reflections
    calculate_bbox(self.predicted_reflections)

  def export_predicted_reflections(self, file_name='predictions.pickle'):
    from dials.model.serialize import dump
    dump.reflections(self.predicted_reflections, file_name)

  def apply_symmetry(self, crystal_model, target_symmetry):
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
      if target_unit_cell is None:
        target_unit_cell = unit_cell
      symm_target_sg = crystal.symmetry(
        unit_cell=target_unit_cell,
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
    from rstbx.symmetry.constraints import parameter_reduction
    s = parameter_reduction.symmetrize_reduce_enlarge(target_space_group)
    s.set_orientation(crystal_model.get_A())
    s.symmetrize()
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
      self, crystal_model, tolerance=0.3):
    print "Candidate crystal model:"
    print crystal_model

    n_rejects = 0

    miller_indices = flex.miller_index()
    self.indexed_reflections = flex.size_t()

    A = crystal_model.get_A()
    A_inv = A.inverse()

    for i_ref in self.reflections_in_scan_range:
      i_rlp = flex.first_index(self.reflections_in_scan_range, i_ref)
      rlp = self.reciprocal_space_points[i_rlp]
      rlp = matrix.col(rlp)
      spot_resolution = 1/rlp.length()
      if spot_resolution < self.d_min:
        continue
      refl = self.reflections[i_ref]
      hkl_float = A_inv * rlp
      hkl_int = [int(round(h)) for h in hkl_float]
      max_difference = max([abs(hkl_float[i] - hkl_int[i]) for i in range(3)])
      if max_difference> tolerance:
        n_rejects += 1
        continue
      miller_indices.append(hkl_int)
      refl.miller_index = hkl_int
      self.indexed_reflections.append(i_ref)

    print "%i reflections indexed successfully (%i rejects)" %(
      self.indexed_reflections.size(), n_rejects)

  def refine(self, crystal_model):
    from dials.algorithms.spot_prediction import ray_intersection
    reflections_for_refinement = ray_intersection(
      self.detector, self.reflections.select(self.indexed_reflections))

    params = self.params.refinement
    from dials.algorithms.refinement import RefinerFactory
    refine = RefinerFactory.from_parameters(self.params, params.verbosity)
    refine.prepare(self.sweep, crystal_model, reflections_for_refinement)
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

  def export_as_json(self, compact=False):
    from dials.model.serialize.dump import crystal as dump_crystal
    from dxtbx.serialize import dump
    with open('crystal.json', 'wb') as f:
      dump_crystal(self.candidate_crystal_models[0], f, compact=compact)
    with open('sweep.json', 'wb') as f:
      dump.imageset(self.sweep, f, compact=compact)

  def export_reflections(self, file_name="indexed.pickle", indexed_only=False):
    reflections = self.reflections
    if indexed_only:
      reflections = reflections.select(self.indexed_reflections)
    with open(file_name, 'wb') as f:
      pickle.dump(reflections, f)

  def export_xds_files(self):
    from dxtbx.serialize import xds
    crystal_model = self.candidate_crystal_models[0]
    A = crystal_model.get_A()
    A_inv = A.inverse()
    real_space_a = A_inv.elems[:3]
    real_space_b = A_inv.elems[3:6]
    real_space_c = A_inv.elems[6:9]
    to_xds = xds.to_xds(self.sweep)
    with open('XDS.INP', 'wb') as f:
      to_xds.XDS_INP(out=f, job_card="XYCORR INIT DEFPIX INTEGRATE CORRECT")
    with open('XPARM.XDS', 'wb') as f:
      to_xds.xparm_xds(
        real_space_a, real_space_b, real_space_c,
        crystal_model.get_space_group().type().number(),
        out=f)

def run(args):
  import time
  from libtbx.phil import command_line
  from dxtbx.imageset import ImageSetFactory

  args = sys.argv[1:]
  cmd_line = command_line.argument_interpreter(master_params=master_phil_scope)
  working_phil, args = cmd_line.process_and_fetch(
      args=args, custom_processor="collect_remaining")
  working_phil.show()

  reflections_filename = args[0]
  sweep_filenames = args[1:]

  sweeps = ImageSetFactory.new(sweep_filenames)
  assert len(sweeps) == 1
  sweep = sweeps[0]
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
