from __future__ import division

import iotbx.phil

master_phil_scope = iotbx.phil.parse("""
d_min = None
  .type = float(value_min=0)
""")

master_params = master_phil_scope.fetch().extract()


def run(args):
  from dials.util.command_line import Importer
  from libtbx.phil import command_line

  importer = Importer(args)
  if len(importer.imagesets) == 0:
    print "Must supply an imageset/sweep!"
    return
  elif len(importer.imagesets) > 1:
    raise RuntimeError("Only one imageset can be processed at a time")
  sweep = importer.imagesets[0]
  assert len(importer.reflections) == 1
  reflections = importer.reflections[0]
  assert len(reflections) > 0
  crystals = importer.crystals
  assert len(importer.crystals) > 0

  args = importer.unhandled_arguments

  cmd_line = command_line.argument_interpreter(master_params=master_phil_scope)
  working_phil = cmd_line.process_and_fetch(args=args)
  working_phil.show()
  params = working_phil.extract()

  from dials.algorithms.indexing.indexer import Indexer

  spots_mm = Indexer._map_spots_pixel_to_mm_rad(
    reflections,
    detector=sweep.get_detector(),
    scan=sweep.get_scan())
  reciprocal_space_points = Indexer._map_centroids_to_reciprocal_space(
    spots_mm=spots_mm,
    detector=sweep.get_detector(),
    beam=sweep.get_beam(),
    goniometer=sweep.get_goniometer())

  from scitbx.array_family import flex
  spots_mm.set_crystal(flex.int(spots_mm.size(), -1))

  from dials.scratch.rjg.index_3D_FFT_simple import index_reflections
  index_reflections(
    spots_mm,
    reciprocal_space_points,
    crystals,
    d_min=params.d_min,
    tolerance=0.3,
    verbosity=1)

  import cPickle as pickle
  pickle.dump(spots_mm, open('indexed.pickle', 'wb'))


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
