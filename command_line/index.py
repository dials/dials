from __future__ import division
from libtbx.phil import command_line
from dials.util.command_line import Importer
from dials.algorithms.indexing.indexer2 import master_phil_scope


def run(args):

  importer = Importer(args)
  if len(importer.imagesets) == 0:
    print "No sweep object could be constructed"
    return
  elif len(importer.imagesets) > 1:
    raise RuntimeError("Only one imageset can be processed at a time")
  sweeps = importer.imagesets
  assert(len(importer.reflections) == 1)
  reflections = importer.reflections[0]
  assert len(reflections) > 0
  args = importer.unhandled_arguments

  sweep = sweeps[0]
  cmd_line = command_line.argument_interpreter(master_params=master_phil_scope)
  working_phil = cmd_line.process_and_fetch(args=args)
  working_phil.show()

  gonio = sweep.get_goniometer()
  detector = sweep.get_detector()
  scan = sweep.get_scan()
  beam = sweep.get_beam()
  print detector
  print scan
  print gonio
  print beam

  params = working_phil.extract()
  if params.method == "3d_fft":
    from dials.algorithms.indexing.fft3d import indexer_fft3d as indexer
  elif params.method == "1d_fft":
    from dials.algorithms.indexing.fft1d import indexer_fft1d as indexer
  elif params.method == "real_space_grid_search":
    from dials.algorithms.indexing.real_space_grid_search \
         import indexer_real_space_grid_search as indexer
  idxr = indexer(reflections, sweep, params=params)
  idxr.index()
  return


if __name__ == '__main__':
  import sys
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run(sys.argv[1:])
