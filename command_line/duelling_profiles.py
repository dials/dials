from __future__ import division

# LIBTBX_SET_DISPATCHER_NAME dev.dials.duelling_profiles

import iotbx.phil

phil_scope = iotbx.phil.parse("""\
show_all_reflection_data = False
  .type = bool
  .help = "Whether or not to print individual reflections"
""", process_includes=True)

help_message = '''

Examples::

  dev.dials.duelling_profiles experiments.json reflections.pickle

'''

def main(reflections, experiments):
  print len(reflections)
  print experiments.crystal.num_scan_points

def run(args):
  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  from dials.util.options import flatten_datablocks
  from dials.util.options import flatten_reflections
  import libtbx.load_env

  usage = "%s [options] datablock.json experiments.json" % (
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_experiments=True,
    read_reflections=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=False)
  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)

  if len(experiments) != 1 or len(reflections) != 1:
    parser.print_help()
    exit()

  main(reflections[0], experiments[0])

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
