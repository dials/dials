from __future__ import division
from dials.util.export_mtz import export_mtz
from libtbx.phil import parse

help_message = '''

This program is used to export the results of dials processing as an
unmerged mtz file, ready for input to downstream programs such as Pointless
and Aimless. The required input is an experiments.json file and an
integrated.pickle file. Optionally the name of the output mtz file can be
specified.

Examples::

  dials.export_mtz experiments.json integrated.pickle

  dials.export_mtz experiments.json integrated.pickle hklout=integrated.mtz

'''

phil_scope = parse('''
  hklout = hklout.mtz
    .type = str
    .help = "The output MTZ file"
  ignore_panels = False
    .type = bool
    .help = "Ignore multiple panels / detectors in output"
  include_partials = False
    .type = bool
    .help = "Include partial reflections (scaled) in output"
  keep_partials = False
    .type = bool
    .help = "Keep low partiality reflections"
  min_isigi = -5
    .type = float
    .help = "Exclude reflections with unfeasible values of I/Sig(I)"

  output {
   
    log = dials.export_mtz.log
      .type = str
    debug_log = dials.export_mtz.debug.log
      .type = str

  }
''')

def run(args):

  import libtbx.load_env
  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  from dials.util.options import flatten_reflections
  from dials.util import log
  from logging import info
  info_handle = log.info_handle()

  usage = '%s integrated.pickle experiments.json [options]' % (
              libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage = usage,
    read_experiments=True,
    read_reflections=True,
    check_format=False,
    phil=phil_scope,
    epilog=help_message)
  params, options = parser.parse_args(show_diff_phil=False)

  # Configure the logging
  log.config(
    info=params.output.log,
    debug=params.output.debug_log)

  from dials.util.version import dials_version
  info(dials_version())

  # Log the diff phil
  diff_phil = parser.diff_phil.as_str()
  if diff_phil is not '':
    info('The following parameters have been modified:\n')
    info(diff_phil)

  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)
  if len(reflections) == 0 or len(experiments) == 0:
    parser.print_help()
    return

  integrated_data = reflections[0]
  experiment_list = experiments
  m = export_mtz(integrated_data, experiment_list, params.hklout,
                 ignore_panels=params.ignore_panels,
                 include_partials=params.include_partials,
                 keep_partials=params.keep_partials,
                 min_isigi=params.min_isigi)
  m.show_summary()


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
