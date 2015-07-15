from __future__ import division
from libtbx.phil import parse

help_message = '''

This program is used to export the results of dials processing as an nxmx file.
The required input is an experiments.json file and an integrated.pickle file.
Optionally the name of the output mtz file can be specified.

Examples::

  dials.export_nxmx experiments.json integrated.pickle

  dials.export_nxmx experiments.json integrated.pickle hklout=integrated.nxs

'''

phil_scope = parse('''
  hklout = hklout.nxs
    .type = str
    .help = "The output NXmx file"
''')

def run(args):
  import libtbx.load_env
  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  from dials.util.options import flatten_reflections
  from dials.util.nexus import dump

  usage = '%s integrated.pickle experiments.json [options]' % (
              libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage = usage,
    read_experiments=True,
    read_reflections=True,
    check_format=False,
    phil=phil_scope,
    epilog=help_message)
  params, options = parser.parse_args(show_diff_phil=True)
  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)
  if len(reflections) == 0 or len(experiments) == 0:
    parser.print_help()
    return

  integrated_data = reflections[0]
  dump(
    experiments,
    integrated_data,
    params.hklout)


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
