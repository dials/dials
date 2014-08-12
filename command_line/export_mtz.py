from __future__ import division
from dials.util.export_mtz import export_mtz

def run(args):
  if len(args) not in (2, 3):
    from libtbx.utils import Usage
    import libtbx.load_env
    usage_message = """\
%s integrated.pickle experiments.json [hklout.mtz] [options]\
""" %libtbx.env.dispatcher_name
    raise Usage(usage_message)

  from dials.util.command_line import Importer
  importer = Importer(args, check_format=False)
  assert importer.experiments is not None and len(importer.experiments) > 0
  assert importer.reflections is not None and len(importer.reflections) == 1

  args = importer.unhandled_arguments
  if len(args) == 0:
    hklout = "hklout.mtz"
  else:
    assert len(args) == 1
    hklout = args[0]

  integrated_data = importer.reflections[0]
  experiment_list = importer.experiments
  m = export_mtz(integrated_data, experiment_list, hklout)
  m.show_summary()


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
