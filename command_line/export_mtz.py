from __future__ import division
from dials.util.export_mtz import export_mtz

def run(args):

  import libtbx.load_env
  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  from dials.util.options import flatten_reflections
  from libtbx.phil import parse
  phil_scope = parse('''
    hklout = hklout.mtz
      .type = str
      .help = "The output MTZ file"
    ignore_panels = False
      .type = bool
      .help = "Ignore multiple panels / detectors in output"
  ''')

  usage = '%s integrated.pickle experiments.json [hklout.mtz] [options]' % (
              libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage = usage,
    read_experiments=True,
    read_reflections=True,
    check_format=False,
    phil=phil_scope)
  params, options = parser.parse_args(show_diff_phil=True)
  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)
  if len(reflections) == 0 or len(experiments) == 0:
    parser.print_help()
    return

  integrated_data = reflections[0]
  experiment_list = experiments
  m = export_mtz(integrated_data, experiment_list, params.hklout,
                 params.ignore_panels)
  m.show_summary()


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
