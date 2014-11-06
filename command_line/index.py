from __future__ import division

try:
  # try importing scipy.linalg before any cctbx modules, otherwise we
  # sometimes get a segmentation fault/core dump if it is imported after
  # scipy.linalg is a dependency of sklearn.cluster.DBSCAN
  import scipy.linalg # import dependency
except ImportError, e:
  pass

from libtbx.phil import command_line
from dials.util.options import OptionParser
from dials.util.options import flatten_reflections
from dials.util.options import flatten_datablocks
from dials.util.options import flatten_experiments


def run(args):
  import libtbx.load_env
  from libtbx.utils import Sorry
  usage = """\
%s [options] datablock.json strong.pickle

Parameters:
""" %libtbx.env.dispatcher_name

  import iotbx.phil
  master_phil_scope = iotbx.phil.parse("""\
include scope dials.algorithms.indexing.indexer.master_phil_scope
output {
  experiments = experiments.json
    .type = path
  reflections = indexed.pickle
    .type = path
}
""", process_includes=True)

  parser = OptionParser(
    usage=usage,
    phil=master_phil_scope,
    read_reflections=True,
    read_datablocks=True,
    read_experiments=True,
    check_format=False)

  params, options = parser.parse_args(show_diff_phil=True)
  datablocks = flatten_datablocks(params.input.datablock)
  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)

  if len(datablocks) == 0:
    if len(experiments) > 0:
      imagesets = importer.experiments.imagesets()
    else:
      parser.print_help()
      return
  elif len(datablocks) > 1:
    raise Sorry("Only one DataBlock can be processed at a time")
  else:
    imagesets = datablocks[0].extract_imagesets()
  if len(experiments):
    known_crystal_models = experiments.crystals()
  else:
    known_crystal_models = None
  assert(len(reflections) == 1)
  reflections = reflections[0]

  for imageset in imagesets:
    if (imageset.get_goniometer() is not None and
        imageset.get_scan() is not None and
        imageset.get_scan().get_oscillation()[1] == 0):
      imageset.set_goniometer(None)
      imageset.set_scan(None)

  if known_crystal_models is not None:
    from dials.algorithms.indexing.known_orientation \
         import indexer_known_orientation
    idxr = indexer_known_orientation(
      reflections, imagesets, params, known_crystal_models)
  elif params.indexing.method == "fft3d":
    from dials.algorithms.indexing.fft3d import indexer_fft3d
    idxr = indexer_fft3d(reflections, imagesets, params=params)
  elif params.indexing.method == "fft1d":
    from dials.algorithms.indexing.fft1d import indexer_fft1d
    idxr = indexer_fft1d(reflections, imagesets, params=params)
  elif params.indexing.method == "real_space_grid_search":
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
