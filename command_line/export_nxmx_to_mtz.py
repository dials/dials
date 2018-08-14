from __future__ import absolute_import, division, print_function

from libtbx.phil import parse

# LIBTBX_SET_DISPATCHER_NAME dev.dials.export_nxms_to_mtz

help_message = '''

This program is used to export an NXmx file as an mtz file.
Optionally the name of the output mtz file can be specified.

Examples::

  dials.export_nxmx_to_mtz hklin=integrated.nxs

  dials.export_nxmx_to_mtz hklin=integrated.nxs hklout=integrated.mtz

'''

phil_scope = parse('''
  hklin = None
    .type = str
    .help = "The input NXmx file"
  hklout = hklout.mtz
    .type = str
    .help = "The output mtz file"
  include scope dials.command_line.export.phil_scope
''', process_includes=True)

def run(args):
  from dials.util.nexus import load
  from dials.util.export_mtz import export_mtz
  import libtbx.load_env
  from dials.util.options import OptionParser

  usage = '%s hklin=hklin.nxs hklout=hklout.mtz [options]' % (
              libtbx.env.dispatcher_name)

  parser = OptionParser(usage=usage, phil=phil_scope, epilog=help_message)
  params, _ = parser.parse_args(show_diff_phil=True)
  params.mtz.hklout = params.hklout

  # Load the experiments and reflections from the NXmx file
  experiments, reflections = load(params.hklin)

  # Export the experiments and reflections to the MTZ file
  m = export_mtz(reflections, experiments, params)
  m.show_summary()

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
