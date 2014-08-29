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
from dials.algorithms.indexing.indexer import master_phil_scope


def run(args):
  if len(args) == 0:
    from libtbx.utils import Usage
    import libtbx.load_env
    from cStringIO import StringIO
    usage_message = """\
%s datablock.json strong.pickle [options]

Parameters:
""" %libtbx.env.dispatcher_name
    s = StringIO()
    master_phil_scope.show(out=s)
    usage_message += s.getvalue()
    raise Usage(usage_message)
  importer = Importer(args, check_format=False)
  if importer.datablocks is None or len(importer.datablocks) == 0:
    if importer.experiments is not None and len(importer.experiments):
      imagesets = importer.experiments.imagesets()
    else:
      print "No DataBlock could be constructed"
      return
  elif len(importer.datablocks) > 1:
    raise RuntimeError("Only one DataBlock can be processed at a time")
  else:
    imagesets = importer.datablocks[0].extract_imagesets()
  if importer.experiments is not None and len(importer.experiments):
    known_crystal_models = importer.experiments.crystals()
  else:
    known_crystal_models = None
  assert len(importer.reflections) == 1
  reflections = importer.reflections[0]
  args = importer.unhandled_arguments

  for imageset in imagesets:
    if (imageset.get_goniometer() is not None and
        imageset.get_scan() is not None and
        imageset.get_scan().get_oscillation()[1] == 0):
      imageset.set_goniometer(None)
      imageset.set_scan(None)

  cmd_line = command_line.argument_interpreter(master_params=master_phil_scope)
  working_phil = cmd_line.process_and_fetch(args=args)
  working_phil.show()

  #filenames = imagesets[0].paths()
  #beam = imagesets[0].get_beam()
  #detector = imagesets[0].get_detector()

  #from dxtbx.imageset import ImageSet
  #from dxtbx.imageset import NullReader
  #reader = NullReader(filenames)
  #imgset = ImageSet(reader)
  #imgset.set_beam(beam)
  #imgset.set_detector(detector)
  #print imgset.get_beam()
  #print imgset.get_detector()
  #imagesets = [imgset]

  params = working_phil.extract()
  if known_crystal_models is not None:
    from dials.algorithms.indexing.known_orientation \
         import indexer_known_orientation
    idxr = indexer_known_orientation(
      reflections, imagesets, params, known_crystal_models)
  elif params.method == "fft3d":
    from dials.algorithms.indexing.fft3d import indexer_fft3d
    idxr = indexer_fft3d(reflections, imagesets, params=params)
  elif params.method == "fft1d":
    from dials.algorithms.indexing.fft1d import indexer_fft1d
    idxr = indexer_fft1d(reflections, imagesets, params=params)
  elif params.method == "real_space_grid_search":
    from dials.algorithms.indexing.real_space_grid_search \
         import indexer_real_space_grid_search
    idxr = indexer_real_space_grid_search(reflections, imagesets, params=params)
  refined_experiments = idxr.refined_experiments
  refined_reflections = idxr.refined_reflections
  if len(refined_experiments):
    idxr.export_as_json(refined_experiments,
                        file_name=params.output.experiments)
    idxr.export_reflections(
      refined_reflections, file_name=params.output.reflections)

  return


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
