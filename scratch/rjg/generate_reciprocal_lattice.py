from __future__ import division
import scitbx.matrix
from scitbx import fftpack
from cctbx.array_family import flex
from cctbx import crystal, sgtbx, uctbx, xray
from dials.model.experiment.crystal_model import Crystal
from dxtbx.serialize.load import imageset_from_string
from libtbx.math_utils import iceil
import iotbx.phil
import math
import sys


master_phil_scope = iotbx.phil.parse("""
unit_cell = None
  .type = unit_cell
space_group = None
  .type = space_group
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
random_seed = 42
  .type = int(value_min=0)
rotation_angle = 90
  .type = float(value_min=0)
""")

class gen_lattice_points(object):
  def __init__(self, sweep, params):
    self.params = params
    flex.set_random_seed(params.random_seed)
    unit_cell = params.unit_cell
    assert unit_cell is not None
    sgi = params.space_group
    if sgi is None:
      sgi = sgtbx.space_group_info(symbol="P 1")
    B = scitbx.matrix.sqr(unit_cell.fractionalization_matrix()).transpose()
    U = scitbx.matrix.sqr(flex.random_double_r3_rotation_matrix())
    direct_matrix = (U * B).inverse()
    crystal_model = Crystal(direct_matrix[0:3],
                            direct_matrix[3:6],
                            direct_matrix[6:9],
                            space_group=sgi.group())
    scan = sweep.get_scan()
    angle = self.params.rotation_angle
    scan.set_image_range((1, iceil(angle/scan.get_oscillation()[1])))
    predicted = predict_reflections(sweep, crystal_model)
    beam_vectors = predicted.beam_vector()
    S = beam_vectors - sweep.get_beam().get_s0()
    centroids = S.rotate_around_origin(sweep.get_goniometer().get_rotation_axis(),
                                       -predicted.rotation_angle())
    self.d_min = self.params.reciprocal_space_grid.d_min
    self.gridding = tuple([self.params.reciprocal_space_grid.n_points]*3)
    centroids = centroids.select((1/centroids.norms())>=self.d_min)
    assert len(centroids) > 0
    self.map_to_grid(sweep, centroids)
    self.fft()
    debug_write_reciprocal_lattice_points_as_pdb(centroids)
    self.debug_write_ccp4_map(self.grid_real, "fft.map")

  def map_to_grid(self, sweep, centroids):
    b_iso = 200
    beam = sweep.get_beam()
    wavelength = beam.get_wavelength()
    d_min = self.d_min

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

    for i_pnt, point in enumerate(centroids):
      point = scitbx.matrix.col(point)
      spot_resolution = 1/point.length()
      if spot_resolution < d_min:
        continue

      grid_coordinates = [int(round(point[i]/rlgrid)+n_points/2) for i in range(3)]
      if max(grid_coordinates) >= n_points: continue # this reflection is outside the grid
      if min(grid_coordinates) < 0: continue # this reflection is outside the grid
      T = math.exp(b_iso * point.length()**2 / 4)
      grid[grid_coordinates] = T

    self.reciprocal_space_grid = grid

  def fft(self):
    fft = fftpack.complex_to_complex_3d(self.gridding)
    grid_complex = flex.complex_double(
      reals=self.reciprocal_space_grid,
      imags=flex.double(self.reciprocal_space_grid.size(), 0))
    grid_transformed = fft.forward(grid_complex)
    #self.grid_real = flex.pow2(flex.abs(grid_transformed))
    self.grid_real = flex.abs(grid_transformed)
    #self.grid_real = flex.pow2(flex.real(grid_transformed))
    #self.grid_real = flex.pow2(flex.imag(self.grid_transformed))
    del grid_transformed

  def debug_write_ccp4_map(self, map_data, file_name):
    from iotbx import ccp4_map
    gridding_first = (0,0,0)
    gridding_last = map_data.all()
    labels = ["cctbx.miller.fft_map"]
    ccp4_map.write_ccp4_map(
      file_name=file_name,
      unit_cell=self.fft_cell,
      space_group=sgtbx.space_group("P1"),
      gridding_first=gridding_first,
      gridding_last=gridding_last,
      map_data=map_data,
      labels=flex.std_string(labels))


def debug_write_reciprocal_lattice_points_as_pdb(
    points, file_name='reciprocal_lattice.pdb'):
  from cctbx import crystal, xray
  cs = crystal.symmetry(unit_cell=(1000,1000,1000,90,90,90), space_group="P1")
  xs = xray.structure(crystal_symmetry=cs)
  sel = flex.random_selection(points.size(),
                              min(20000, points.size()))
  rsp = points.select(sel)
  for site in rsp:
    xs.add_scatterer(xray.scatterer("C", site=site))

  xs.sites_mod_short()
  with open(file_name, 'wb') as f:
    print >> f, xs.as_pdb_file()


def predict_reflections(sweep, crystal_model):
  from dials.algorithms.integration import ReflectionPredictor
  predictor = ReflectionPredictor()
  reflections = predictor(sweep, crystal_model)
  return reflections


def run(args):
  from libtbx.phil import command_line
  from dials.util.command_line import Importer

  args = sys.argv[1:]
  importer = Importer(args)
  if len(importer.imagesets) == 0:
    print "No sweep object could be constructed"
    return
  elif len(importer.imagesets) > 1:
    raise RuntimeError("Only one imageset can be processed at a time")
  sweeps = importer.imagesets
  args = importer.unhandled_arguments

  sweep = sweeps[0]
  cmd_line = command_line.argument_interpreter(master_params=master_phil_scope)
  working_phil, args = cmd_line.process_and_fetch(
      args=args, custom_processor="collect_remaining")
  working_phil.show()

  result = gen_lattice_points(sweep, params=working_phil.extract())


if __name__ == '__main__':
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run(sys.argv[1:])
