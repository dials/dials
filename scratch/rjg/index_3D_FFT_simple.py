from __future__ import division
import sys
import cPickle as pickle
import math

import libtbx.phil
from scitbx import matrix
from scitbx import fftpack

from cctbx.array_family import flex
from cctbx import crystal, miller, sgtbx, uctbx, xray

import dxtbx
from dials.model.data import ReflectionList
from dials.model.experiment.crystal_model import Crystal


master_phil_scope = libtbx.phil.parse("""
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
debug = False
  .type = bool
refinement {
  n_macro_cycles = 3
    .type = int(value_min=0)
  fix_detector = False
    .type = bool
    .help = "Whether or not to refine the detector position and orientation."
  fix_beam = False
    .type = bool
    .help = "Whether or not to refine the beam direction."
}
""")

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
  assert length_b > length_a
  angle = vec_a.angle(vec_b, deg=True)
  if angle < angular_tolerance or abs(180-angle) < angular_tolerance:
    n = length_b/length_a
    if abs(round(n) - n) < relative_tolerance:
      return True
  return False


class indexer(object):

  def __init__(self, reflections, goniometer, detector, scan, beam,
               params=None):
    self.reflections = reflections
    self.goniometer = goniometer
    self.detector = detector
    self.scan = scan
    self.beam = beam
    if params is None: params = master_params
    self.params = params

  def index(self):
    n_points = self.params.fft_n_points
    self.gridding = fftpack.adjust_gridding_triple(
      (n_points,n_points,n_points), max_prime=5)
    n_points = self.gridding[0]
    self.map_centroids_to_reciprocal_space_grid()

    self.unit_cell = uctbx.unit_cell([n_points*self.params.d_min/2]*3+[90]*3)
    self.crystal_symmetry = crystal.symmetry(unit_cell=self.unit_cell,
                                             space_group_symbol="P1")

    print "Number of centroids used: %i" %(
      (self.reciprocal_space_grid>0).count(True))
    self.fft()
    self.find_peaks()
    self.find_candidate_basis_vectors()
    self.find_candidate_orientation_matrices()
    self.index_reflections_given_orientation_matix(
      self.candidate_crystal_models[0])
    for i in range(self.params.refinement.n_macro_cycles):
      print "Starting refinement (macro-cycle %i)" %(i+1)
      print
      self.refine(self.candidate_crystal_models[0])
      self.index_reflections_given_orientation_matix(
        self.candidate_crystal_models[0])

    if self.params.debug:
      self.debug_write_reciprocal_lattice_points_as_pdb()
      self.debug_write_ccp4_map(map_data=self.grid_real, file_name="patt.map")
      self.debug_show_candidate_basis_vectors()

  def map_centroids_to_reciprocal_space_grid(self):

    wavelength = self.beam.get_wavelength()
    s0 = matrix.col(self.beam.get_s0())
    rotation_axis = matrix.col(
      self.goniometer.get_rotation_axis())

    self.reciprocal_space_points = []

    n_points = self.gridding[0]
    rlgrid = 2 * wavelength / (self.params.d_min * n_points)
    print "FFT gridding: (%i,%i,%i)" %self.gridding

    grid = flex.double(flex.grid(self.gridding), 0)

    reflections_used_for_indexing = ReflectionList()

    deg_to_radians = math.pi/180

    oscillation_range = self.scan.get_oscillation_range(deg=True)
    phi_zero = oscillation_range[0]
    image_width_rad = (
      oscillation_range[1] - oscillation_range[0]) * deg_to_radians

    for refl in self.reflections:

      frame_number = refl.centroid_position[2]

      if len(self.params.scan_range):
        use_reflection = False
        for scan_range in self.params.scan_range:
          if scan_range is None: continue
          range_start, range_end = scan_range
          if frame_number >= range_start and frame_number < range_end:
            use_reflection = True
            break
        if not use_reflection:
          continue

      s1 = matrix.col(
        self.detector.get_pixel_lab_coord(refl.centroid_position[:2]))
      s1 = s1.normalize()/wavelength
      S = s1 - s0
      spot_resolution = 1/S.length()

      if spot_resolution < self.params.d_min:
        continue

      phi = self.scan.get_angle_from_array_index(frame_number)

      # set reflection properties that might be needed by the dials refinement
      # engine, and convert values from pixels and image number to mm/rads
      from dials.algorithms.centroid import centroid_px_to_mm
      refl.beam_vector = tuple(s1)
      refl.frame_number = refl.centroid_position[2]
      centroid_position, centroid_variance, _ = centroid_px_to_mm(
        self.detector, self.scan,
        refl.centroid_position,
        refl.centroid_variance,
        (1,1,1))
      refl.beam_vector = tuple(s1)
      refl.frame_number = refl.centroid_position[2]
      refl.centroid_position = centroid_position
      refl.centroid_variance = centroid_variance
      refl.rotation_angle = centroid_position[2]

      point = S.rotate_around_origin(rotation_axis, -phi, deg=True)

      #assert round(detector.get_resolution_at_pixel(beam.get_s0(), wavelength,
                                                    #refl.centroid_position[:2]),
                   #4) == \
             #round(spot_resolution, 4)

      grid_coordinates = [int(round(point[i]/rlgrid)+n_points/2) for i in range(3)]
      if max(grid_coordinates) >= n_points: continue # this reflection is outside the grid
      if min(grid_coordinates) < 0: continue # this reflection is outside the grid
      T = math.exp(-self.params.b_iso * S.length()**2 / 4)
      grid[grid_coordinates] = T
      self.reciprocal_space_points.append(tuple(point))
      reflections_used_for_indexing.append(refl)

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
    from cctbx import masks
    flood_fill = masks.flood_fill(grid_real_binary.iround(), self.unit_cell)
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
    length_tolerance = 5 # Angstrom
    # XXX maybe the tolerance should be unitless, i.e. proportional to vector length?
    angle_tolerance = 10 # degrees

    orth = self.unit_cell.orthogonalize
    for v in vectors:
      length = v.length()
      if length < self.params.min_cell or length > self.params.max_cell:
        continue
      matched_group = False
      for group in vector_groups:
        mean_v = group.mean()
        if abs(mean_v.length() - length) < length_tolerance:
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
          model = Crystal(a, b, c)
          uc = model.get_unit_cell()
          params = uc.parameters()
          if uc.volume() > (params[0]*params[1]*params[2]/100):
            # unit cell volume cutoff from labelit 2004 paper
            self.candidate_crystal_models.append(model)

  def index_reflections_given_orientation_matix(
      self, crystal_model, tolerance=0.3):
    print "Candidate crystal model:"
    print crystal_model

    n_rejects = 0

    miller_indices = flex.miller_index()
    self.indexed_reflections = ReflectionList()

    A = crystal_model.get_A()
    A_inv = A.inverse()

    for rlp, refl in zip(self.reciprocal_space_points, self.reflections_used_for_indexing):
      hkl_float = A_inv * matrix.col(rlp)
      hkl_int = [int(round(h)) for h in hkl_float]
      max_difference = max([abs(hkl_float[i] - hkl_int[i]) for i in range(3)])
      if max_difference> tolerance:
        n_rejects += 1
        continue
      miller_indices.append(hkl_int)
      refl.miller_index = hkl_int
      self.indexed_reflections.append(refl)

    print "%i reflections indexed successfully (%i rejects)" %(
      self.indexed_reflections.size(), n_rejects)

  def refine(self, crystal_model):
    from  dials.algorithms.refinement import refine
    from dials.algorithms.spot_prediction import ray_intersection
    indexed_reflections = ray_intersection(
      self.detector, self.indexed_reflections)

    print "Starting crystal model:"
    print crystal_model

    print "Starting detector model:"
    print self.detector

    print "Starting beam model:"
    print self.beam

    refine(self.beam, self.goniometer, crystal_model, self.detector, self.scan,
           indexed_reflections, verbosity=1,
           fix_cell=False,
           fix_beam=self.params.refinement.fix_beam,
           fix_detector=self.params.refinement.fix_detector,
           scan_varying=False)

    print "Refined crystal model:"
    print crystal_model

    print "Refined detector model:"
    print self.detector

    print "Refined beam model:"
    print self.beam

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

def run(args):
  from libtbx.phil import command_line
  import time

  args = sys.argv[1:]
  cmd_line = command_line.argument_interpreter(master_params=master_phil_scope)
  working_phil, args = cmd_line.process_and_fetch(
      args=args, custom_processor="collect_remaining")
  working_phil.show()

  reflections_filename = args[0]
  sweep_filenames = args[1:]

  models = dxtbx.load(sweep_filenames[0])
  gonio = models.get_goniometer()
  detector = models.get_detector()
  scan = models.get_scan()
  # the refinement MUST have the correct image/oscillation ranges set!
  # XXX not sure how safe this is? should load a sweep instead
  scan.set_image_range((1, len(sweep_filenames)))
  beam = models.get_beam()
  print detector
  print scan
  print gonio
  print beam

  t1 = time.time()
  with open(reflections_filename, 'rb') as f:
    reflections = pickle.load(f)
  t2 = time.time()
  print "Time taken loading reflection file: %.3fs" %(t2-t1)

  idxr = indexer(reflections, gonio, detector, scan, beam,
                 params=working_phil.extract())
  idxr.index()
  return


if __name__ == '__main__':
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run(sys.argv[1:])
