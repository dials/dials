import iotbx.phil

master_params = ""
master_phil_scope = iotbx.phil.parse(input_string=master_params)

def run(args):
  import time
  from libtbx.phil import command_line
  from dials.util.command_line import Importer

  args = sys.argv[1:]
  importer = Importer(args)
  if len(importer.imagesets) == 0:
      print "No sweep object could be constructed"
      return
  elif len(importer.imagesets) > 1:
      raise RuntimeError("Only one imageset can be processed at a time")
  sweeps = importer.imagesets
  reflections = importer.reflections
  assert len(reflections) > 0
  args = importer.unhandled_arguments

  sweep = sweeps[0]
  cmd_line = command_line.argument_interpreter(master_params=master_phil_scope)
  working_phil, args = cmd_line.process_and_fetch(
      args=args, custom_processor="collect_remaining")
  working_phil.show()

  from rstbx.phil.phil_preferences import libtbx_defs
  import iotbx.phil
  hardcoded_phil = iotbx.phil.parse(input_string=libtbx_defs).extract()

  gonio = sweep.get_goniometer()
  detector = sweep.get_detector()
  scan = sweep.get_scan()
  beam = sweep.get_beam()
  print detector
  print scan
  print gonio
  print beam

  from dials.algorithms.indexing import indexer

  #strategies = None
  #idxr = indexer.Indexer(strategies, parameters=working_phil.extract())
  #idxr(reflections, detector, beam, goniometer=gonio, scan=scan)

  # first exercise the better experimental model discovery

  detector, beam = indexer.discover_better_experimental_model(
    reflections, detector, beam,
    goniometer=gonio, scan=scan, params=hardcoded_phil)

  print detector
  print beam


if __name__ == '__main__':
  import sys
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run(sys.argv[1:])
