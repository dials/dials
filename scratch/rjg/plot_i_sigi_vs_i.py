from __future__ import division

import iotbx.phil
from dials.util.options import OptionParser

help_message = '''
'''

phil_scope = iotbx.phil.parse("""
""", process_includes=True)


def run(args):
  import libtbx.load_env
  usage = "%s [options]" %libtbx.env.dispatcher_name

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    check_format=False,
    epilog=help_message)

  params, options, args = parser.parse_args(show_diff_phil=True,
                                            return_unhandled=True)

  assert len(args) == 1
  from iotbx.reflection_file_reader import any_reflection_file

  intensities = None

  f = args[0]

  arrays = any_reflection_file(f).as_miller_arrays(merge_equivalents=False)
  for ma in arrays:
    print ma.info().labels
    if ma.info().labels == ['I', 'SIGI']:
      intensities = ma
    elif ma.info().labels == ['IMEAN', 'SIGIMEAN']:
      intensities = ma
    elif ma.info().labels == ['I(+)', 'SIGI(+)', 'I(-)', 'SIGI(-)']:
      intensities = ma

  assert intensities is not None

  i_sigi = intensities.data()/intensities.sigmas()

  # set backend before importing pyplot
  import matplotlib
  matplotlib.use('Agg')

  from matplotlib import pyplot
  pyplot.scatter(intensities.data(), i_sigi, marker='+', s=2, alpha=0.5, c='black')
  pyplot.gca().set_xscale('log',basex=10)
  xlim = pyplot.xlim()
  pyplot.xlim(1, xlim[1])
  pyplot.savefig('i_sigi_vs_i.png')
  pyplot.clf()

  return


if __name__ == '__main__':
  import sys
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run(sys.argv[1:])
