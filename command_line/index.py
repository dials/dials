from __future__ import division

try:
  # try importing scipy.linalg before any cctbx modules, otherwise we
  # sometimes get a segmentation fault/core dump if it is imported after
  # scipy.linalg is a dependency of sklearn.cluster.DBSCAN
  import scipy.linalg # import dependency
except ImportError, e:
  pass

from libtbx.phil import command_line
from dials.util.command_line import Importer
from dials.algorithms.indexing.indexer2 import master_phil_scope


def run(args):

  importer = Importer(args, check_format=False)
  if len(importer.datablocks) == 0:
    print "No DataBlock could be constructed"
    return
  elif len(importer.datablocks) > 1:
    raise RuntimeError("Only one DataBlock can be processed at a time")
  imagesets = importer.datablocks[0].extract_imagesets()
  if len(imagesets) > 1:
    raise RuntimeError("Only one imageset can be processed at a time")
  imageset = imagesets[0]
  assert len(importer.reflections) == 1
  reflections = importer.reflections[0]
  args = importer.unhandled_arguments

  cmd_line = command_line.argument_interpreter(master_params=master_phil_scope)
  working_phil = cmd_line.process_and_fetch(args=args)
  working_phil.show()

  gonio = imageset.get_goniometer()
  detector = imageset.get_detector()
  scan = imageset.get_scan()
  beam = imageset.get_beam()
  print detector
  print scan
  print gonio
  print beam

  params = working_phil.extract()
  if params.method == "fft3d":
    from dials.algorithms.indexing.fft3d import indexer_fft3d as indexer
  elif params.method == "fft1d":
    from dials.algorithms.indexing.fft1d import indexer_fft1d as indexer
  elif params.method == "real_space_grid_search":
    from dials.algorithms.indexing.real_space_grid_search \
         import indexer_real_space_grid_search as indexer
  idxr = indexer(reflections, imageset, params=params)
  idxr.index()
  refined_experiments = idxr.refined_experiments
  refined_reflections = idxr.refined_reflections
  if len(refined_experiments):
    idxr.export_as_json(refined_experiments,
                        file_name=params.output.experiments_filename)
    idxr.export_reflections(
      refined_reflections, file_name=params.output.reflections_filename)

  return


if __name__ == '__main__':
  import sys
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run(sys.argv[1:])
