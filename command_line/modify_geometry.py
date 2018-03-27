from __future__ import absolute_import, division
from __future__ import print_function
import libtbx.phil

help_message = '''

'''

phil_scope= libtbx.phil.parse("""
include scope dials.util.options.geometry_phil_scope
output {
  datablock = modified_datablock.json
    .type = path
  experiments = modified_experiments.json
    .type = path
}
""", process_includes=True)


def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_datablocks
  from dials.util.options import flatten_experiments
  import libtbx.load_env

  usage = "%s [options] datablock.json | experiments.json" %(
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_datablocks=True,
    read_experiments=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)
  experiments = flatten_experiments(params.input.experiments)
  datablocks = flatten_datablocks(params.input.datablock)

  if len(experiments) == 0 and len(datablocks) == 0:
    parser.print_help()
    exit(0)

  from dials.command_line.dials_import import ManualGeometryUpdater
  update_geometry = ManualGeometryUpdater(params)

  if len(experiments):
    imagesets = experiments.imagesets()

  elif len(datablocks):

    assert len(datablocks) == 1
    imagesets = datablocks[0].extract_imagesets()

  for imageset in imagesets:
    imageset_new = update_geometry(imageset)
    imageset.set_detector(imageset_new.get_detector())
    imageset.set_beam(imageset_new.get_beam())
    imageset.set_goniometer(imageset_new.get_goniometer())
    imageset.set_scan(imageset_new.get_scan())

  from dxtbx.serialize import dump
  if len(experiments):
    print("Saving modified experiments to %s" %params.output.experiments)
    dump.experiment_list(experiments, params.output.experiments)
  elif len(datablocks):
    print("Saving modified datablock to %s" %params.output.datablock)
    dump.datablock(datablocks, params.output.datablock)

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
