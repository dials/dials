# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1
# LIBTBX_SET_DISPATCHER_NAME dev.dials.check_strategy

from __future__ import absolute_import, division, print_function

import libtbx.phil
from scitbx.array_family import flex

help_message = '''

'''

phil_scope= libtbx.phil.parse("""
""")


def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  from dials.util.options import flatten_reflections
  import libtbx.load_env

  usage = "%s [options] datablock.json" %(
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_experiments=True,
    read_reflections=True,
    check_format=True,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)
  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)

  if len(experiments) == 0 or len(reflections) == 0:
    parser.print_help()
    exit(0)

  from dials.algorithms.shadowing.filter import filter_shadowed_reflections

  imagesets = experiments.imagesets()
  reflections = reflections[0]
  shadowed = filter_shadowed_reflections(experiments, reflections)

  print("# shadowed reflections: %i/%i (%.2f%%)" %(
    shadowed.count(True), shadowed.size(),
    shadowed.count(True)/shadowed.size() * 100))

  expt = experiments[0]
  x,y,z = reflections['xyzcal.px'].parts()
  z_ = z * expt.scan.get_oscillation()[1]
  zmin, zmax = expt.scan.get_oscillation_range()

  hist_scan_angle = flex.histogram(z_.select(shadowed), n_slots=int(zmax-zmin))
  #hist_scan_angle.show()

  uc = experiments[0].crystal.get_unit_cell()
  d_spacings = uc.d(reflections['miller_index'])
  ds2 = uc.d_star_sq(reflections['miller_index'])

  hist_res = flex.histogram(ds2.select(shadowed), flex.min(ds2), flex.max(ds2), n_slots=20, )
  #hist_res.show()

  import matplotlib
  matplotlib.use('Agg')
  from matplotlib import pyplot as plt

  plt.hist2d(z_.select(shadowed).as_numpy_array(),
             ds2.select(shadowed).as_numpy_array(), bins=(40,40),
             range=((flex.min(z_), flex.max(z_)),(flex.min(ds2), flex.max(ds2))))
  yticks_dsq = flex.double(plt.yticks()[0])
  from cctbx import uctbx
  yticks_d = uctbx.d_star_sq_as_d(yticks_dsq)
  plt.axes().set_yticklabels(['%.2f' %y for y in yticks_d])
  plt.xlabel('Scan angle (degrees)')
  plt.ylabel('Resolution (A^-1)')
  cbar = plt.colorbar()
  cbar.set_label('# shadowed reflections')
  plt.savefig('n_shadowed_hist2d.png')
  plt.clf()

  plt.scatter(hist_scan_angle.slot_centers().as_numpy_array(), hist_scan_angle.slots().as_numpy_array())
  plt.xlabel('Scan angle (degrees)')
  plt.ylabel('# shadowed reflections')
  plt.savefig("n_shadowed_vs_scan_angle.png")
  plt.clf()

  plt.scatter(hist_res.slot_centers().as_numpy_array(), hist_res.slots().as_numpy_array())
  plt.xlabel('d_star_sq')
  plt.savefig("n_shadowed_vs_resolution.png")
  plt.clf()

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
