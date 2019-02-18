# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

# DIALS_ENABLE_COMMAND_LINE_COMPLETION
from __future__ import absolute_import, division

try:
  import matplotlib
except ImportError:
  exit() # To pass through the "make" step, for graphics-free HPC build

# Offline backend
matplotlib.use("Agg")

import logging
logger = logging.getLogger('dials.command_line.rl_png')

import libtbx.phil
from scitbx import matrix
from scitbx.array_family import flex

from dials.command_line.reciprocal_lattice_viewer import render_3d
from dials.command_line.reciprocal_lattice_viewer import help_message


phil_scope= libtbx.phil.parse("""
include scope dials.command_line.reciprocal_lattice_viewer.phil_scope
basis_vector_search {
  n_solutions = 3
    .type = int
}
plot {
  size_inches = 10,10
    .type = floats(size=2, value_min=0)
}
""", process_includes=True)

def settings():
  return phil_scope.fetch().extract()


class ReciprocalLatticePng(render_3d):

  def __init__(self, settings=None):
    render_3d.__init__(self)
    if settings is not None:
      self.settings = settings
    else:
      self.settings = settings()
    self.viewer = PngScene(settings=self.settings)

class PngScene(object):

  def __init__(self, settings):
    self.settings = settings
    self.rotation_axis = None
    self.beam_vector = None
    self.points = None
    self.colors = None
    self.palette = None

  def set_rotation_axis(self, axis):
    self.rotation_axis = axis

  def set_beam_vector(self, beam):
    self.beam_vector = beam

  def set_points(self, points):
    self.points = points

  def set_colors(self, colors):
    import math
    # convert whites to black (background is white)
    colors.set_selected((colors.norms() == math.sqrt(3)), (0,0,0))
    self.colors = colors

  def set_palette(self, palette):
    self.palette = palette

  def project_2d(self, n):
    d = self.points.dot(n.elems)
    p = d * flex.vec3_double(len(d), n.elems)

    points2d = self.points - p

    x = matrix.col((1,0,0))
    if x.angle(n) == 0 or x.angle(-n) == 0:
      x = matrix.col((0,1,0))

    x = (x - x.dot(n) * n).normalize()

    y = x.cross(n)
    #assert y.angle(x, deg=True) == 90
    #assert y.angle(matrix.col(n), deg=True) == 90

    px2d = points2d.dot(x)
    py2d = points2d.dot(y)

    return px2d, py2d

  def plot(self, filename, n=(1,0,0)):
    from matplotlib import pyplot

    n = matrix.col(n).normalize()
    x, y = self.project_2d(n)
    fig = pyplot.figure(figsize=self.settings.plot.size_inches)
    pyplot.scatter(x.as_numpy_array(), y.as_numpy_array(),
                   marker='+', s=self.settings.marker_size, c=list(self.colors))
    pyplot.title('Plane normal: (%.2g, %.2g, %.2g)' %(n.elems))
    fig.savefig(filename)
    pyplot.close()

def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  from dials.util.options import flatten_reflections
  from dials.util import log
  import libtbx.load_env

  usage = "%s [options] experiments.json reflections.pickle" %(
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_experiments=True,
    read_reflections=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args()
  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)

  if len(experiments) == 0 or len(reflections) == 0:
    parser.print_help()
    exit(0)

  # Configure the logging
  log.config(info='dials.rl_png.log')

  # Log the diff phil
  diff_phil = parser.diff_phil.as_str()
  if diff_phil is not '':
    logger.info('The following parameters have been modified:\n')
    logger.info(diff_phil)

  reflections = reflections[0]

  imagesets = experiments.imagesets()

  f = ReciprocalLatticePng(settings=params)
  f.load_models(imagesets, reflections, None)

  imageset = imagesets[0]
  rotation_axis = matrix.col(imageset.get_goniometer().get_rotation_axis())
  s0 = matrix.col(imageset.get_beam().get_s0())

  e1 = rotation_axis.normalize()
  e2 = s0.normalize()
  e3 = e1.cross(e2).normalize()
  #print e1
  #print e2
  #print e3

  f.viewer.plot('rl_rotation_axis.png', n=e1.elems)
  f.viewer.plot('rl_beam_vector', n=e2.elems)
  f.viewer.plot('rl_e3.png', n=e3.elems)

  n_solutions = params.basis_vector_search.n_solutions

  if len(experiments):
    for i, c in enumerate(experiments.crystals()):
      A = matrix.sqr(c.get_A())
      astar = A[:3]
      bstar = A[3:6]
      cstar = A[6:9]

      direct_matrix = A.inverse()
      a = direct_matrix[:3]
      b = direct_matrix[3:6]
      c = direct_matrix[6:9]

      prefix = ''
      if len(experiments.crystals()) > 1:
        prefix = '%i_' %(i+1)

      f.viewer.plot('rl_%sa.png' %prefix, n=a)
      f.viewer.plot('rl_%sb.png' %prefix, n=b)
      f.viewer.plot('rl_%sc.png' %prefix, n=c)

  elif n_solutions:
    from dials.command_line.search_beam_position \
         import run_dps, dps_phil_scope

    hardcoded_phil = dps_phil_scope.extract()
    hardcoded_phil.d_min = params.d_min

    imageset = imagesets[0]

    if 'imageset_id' not in reflections:
      reflections['imageset_id'] = reflections['id']

    reflections.centroid_px_to_mm(
      imageset.get_detector(), scan=imageset.get_scan())

    reflections.map_centroids_to_reciprocal_space(
      detector=imageset.get_detector(), beam=imageset.get_beam(),
      goniometer=imageset.get_goniometer())

    if params.d_min is not None:
      d_spacings = 1/reflections['rlp'].norms()
      sel = d_spacings > params.d_min
      reflections = reflections.select(sel)

    # derive a max_cell from mm spots

    from dials.algorithms.indexing.indexer import find_max_cell
    max_cell = find_max_cell(reflections, max_cell_multiplier=1.3,
                             step_size=45).max_cell

    result = run_dps((imageset, reflections, max_cell, hardcoded_phil))
    solutions = [matrix.col(v) for v in result['solutions']]
    for i in range(min(n_solutions, len(solutions))):
      v = solutions[i]
      #if i > 0:
        #for v1 in solutions[:i-1]:
          #angle = v.angle(v1, deg=True)
          #print angle
      f.viewer.plot('rl_solution_%s.png' %(i+1), n=v.elems)


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
