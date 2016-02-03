# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

# DIALS_ENABLE_COMMAND_LINE_COMPLETION
from __future__ import division

import matplotlib

# Offline backend
matplotlib.use("Agg")

from logging import info

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
""", process_includes=True)

def settings () :
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

  def plot(self, filename, n=(1,0,0)):

    n = matrix.col(n).normalize()
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

    from matplotlib import pyplot
    fig = pyplot.figure(figsize=(10,10))
    pyplot.scatter(
      px2d, py2d, marker='+', s=self.settings.marker_size, c=self.colors)
    pyplot.title('Plane normal: (%.2g, %.2g, %.2g)' %(n.elems))
    fig.savefig(filename)
    pyplot.close()

def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_datablocks
  from dials.util.options import flatten_experiments
  from dials.util.options import flatten_reflections
  from dials.util import log
  import libtbx.load_env

  usage = "%s [options] datablock.json reflections.pickle" %(
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_datablocks=True,
    read_experiments=True,
    read_reflections=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args()
  datablocks = flatten_datablocks(params.input.datablock)
  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)

  if (len(datablocks) == 0 and len(experiments) == 0) or len(reflections) == 0:
    parser.print_help()
    exit(0)

  # Configure the logging
  log.config(info='dials.rl_png.log')

  # Log the diff phil
  diff_phil = parser.diff_phil.as_str()
  if diff_phil is not '':
    info('The following parameters have been modified:\n')
    info(diff_phil)

  reflections = reflections[0]

  if len(datablocks) == 0 and len(experiments) > 0:
    imagesets = experiments.imagesets()
  else:
    imagesets = []
    for datablock in datablocks:
      imagesets.extend(datablock.extract_imagesets())

  f = ReciprocalLatticePng(settings=params)
  f.load_models(imagesets, reflections)

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
      A = c.get_A()
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
    from dials.command_line.discover_better_experimental_model \
         import run_dps, dps_phil_scope

    hardcoded_phil = dps_phil_scope.extract()
    hardcoded_phil.d_min = params.d_min
    result = run_dps((imagesets[0], reflections, hardcoded_phil))
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
