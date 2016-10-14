# LIBTBX_SET_DISPATCHER_NAME dev.dials.best_out

from __future__ import division

import iotbx.phil

help_message = '''

'''

phil_scope = iotbx.phil.parse("""\
n_bins = 100
  .type = int
frames = None
  .type = int
  .multiple = True
""", process_includes=True)

def main():
  import sys
  run(sys.argv[1:])

def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  from dials.util.options import flatten_reflections
  import libtbx.load_env
  from dials.util import best

  usage = "%s [options] experiments.json integrated.pickle" % (
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_experiments=True,
    read_reflections=True,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)
  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)

  if len(experiments) == 0 or len(reflections) == 0:
    parser.print_help()
    exit()

  #assert(len(experiments) == 1)

  experiment = experiments[0]
  reflections = reflections[0]
  imageset = experiment.imageset

  best.write_background_file('bestfile.dat', imageset, n_bins=params.n_bins)

  best.write_integrated_hkl('bestfile', reflections)

  best.write_par_file('bestfile.par', experiment)

if __name__ == '__main__':
  main()

