# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

# DIALS_ENABLE_COMMAND_LINE_COMPLETION
from __future__ import division
from scitbx.array_family import flex

import matplotlib

# Offline backend
matplotlib.use("Agg")

from dials.command_line.reciprocal_lattice_viewer import render_3d
from dials.command_line.reciprocal_lattice_viewer import help_message, phil_scope
from dials.command_line.reciprocal_lattice_viewer import settings

class ReciprocalLatticePng(render_3d):

  def __init__(self, settings=None):
    render_3d.__init__(self)
    self.viewer = PngScene()
    if settings is not None:
      self.settings = settings
    else:
      self.settings = settings()

class PngScene(object):

  def __init__(self):
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
    colors.set_selected((colors.norms() == 1), (0,0,0))
    self.colors = colors

  def plot(self, filename, n=(1,0,0)):

    from scitbx import matrix
    n = matrix.col(n).normalize()
    d = self.points.dot(n.elems)
    p = d * flex.vec3_double(len(d), n.elems)

    points2d = self.points - p

    x = matrix.col((1,0,0))
    if x.angle(n) == 0:
      x = matrix.col((0,1,0))

    x = (x - x.dot(n) * n).normalize()

    y = x.cross(n)
    #assert y.angle(x, deg=True) == 90
    #assert y.angle(matrix.col(n), deg=True) == 90

    px2d = points2d.dot(x)
    py2d = points2d.dot(y)

    from matplotlib import pyplot
    fig = pyplot.figure(figsize=(10,10))
    pyplot.scatter(px2d, py2d, marker='+', s=5, c=self.colors)
    fig.savefig(filename)
    pyplot.close()

def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_datablocks
  from dials.util.options import flatten_experiments
  from dials.util.options import flatten_reflections
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

  params, options = parser.parse_args(show_diff_phil=True)
  datablocks = flatten_datablocks(params.input.datablock)
  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)

  if (len(datablocks) == 0 and len(experiments) == 0) or len(reflections) == 0:
    parser.print_help()
    exit(0)

  reflections = reflections[0]

  if len(datablocks) == 0 and len(experiments) > 0:
    imagesets = experiments.imagesets()
  else:
    imagesets = []
    for datablock in datablocks:
      imagesets.extend(datablock.extract_imagesets())

  f = ReciprocalLatticePng(settings=params)
  f.load_models(imagesets, reflections)
  f.viewer.plot('rl_100.png', n=(1,0,0))
  f.viewer.plot('rl_010.png', n=(0,1,0))
  f.viewer.plot('rl_001.png', n=(0,0,1))

  if len(experiments):
    c = experiments.crystals()[0]
    A = c.get_A()
    astar = A[:3]
    bstar = A[3:6]
    cstar = A[6:9]

    direct_matrix = A.inverse()
    a = direct_matrix[:3]
    b = direct_matrix[3:6]
    c = direct_matrix[6:9]

    f.viewer.plot('rl_a.png', n=a)
    f.viewer.plot('rl_b.png', n=b)
    f.viewer.plot('rl_c.png', n=c)


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
