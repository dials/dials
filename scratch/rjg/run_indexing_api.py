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
  if len(sys.argv) == 1:
    raise RuntimeError, '%s sweep.json reflections.pickle' % sys.argv[0]
  run(sys.argv[1:])

howto = '''
dials.import ~/data/12287/12287_1_E1_0*img
Importing data from the following sources:
 - Sweep from image file data
 - Crystal from nowhere!
 - Parameters from system parameters
Saved sweep to sweep.json
Saved parameters to param.phil

dials.spotfinder scan_range=0,5 scan_range=55,60  ~/data/12287/12287_1_E1_0*img -o refl.pkl
Configuring spot finder from input parameters
Finding strong spots

Finding spots in image 0 to 5...
Extracted 147376 strong pixels.......................................4.89s
Extracted 7487 spots from pixels.....................................0.15s

Finding spots in image 55 to 60...
Extracted 143364 strong pixels.......................................4.09s
Extracted 6515 spots from pixels.....................................0.12s
Calculated 14002 spot centroids......................................0.33s
Calculated 14002 spot intensities....................................0.07s
Filtered 7868 spots by number of pixels..............................0.11s
Filtered 7432 spots by peak-centroid distance........................0.10s
Saved 7432 reflections to refl.pkl...................................0.06s

'''
