from __future__ import division

from libtbx.phil import command_line
import iotbx.phil
from dials.util.command_line import Importer

import libtbx.load_env
dials_path = libtbx.env.dist_path('dials')

master_phil_scope = iotbx.phil.parse("""
include file %s/data/refinement.phil
verbosity = 0
  .type = int(value_min=0)
nproc = Auto
  .type = int(value_min=1)
""" %dials_path, process_includes=True)

master_params = master_phil_scope.fetch().extract()


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
  crystals = importer.crystals
  assert len(crystals) == 1
  args = importer.unhandled_arguments

  sweep = sweeps[0]
  cmd_line = command_line.argument_interpreter(master_params=master_phil_scope)
  working_phil = cmd_line.process_and_fetch(args=args)
  working_phil.show()
  params = working_phil.extract()

  goniometer = sweep.get_goniometer()
  detector = sweep.get_detector()
  scan = sweep.get_scan()
  beam = sweep.get_beam()

  from dials.algorithms.indexing.symmetry \
       import refined_settings_factory_from_refined_triclinic

  Lfat = refined_settings_factory_from_refined_triclinic(
    params, scan, goniometer, beam, detector, crystals[0], reflections,
    nproc=params.nproc, refiner_verbosity=params.verbosity)
  Lfat.labelit_printout()
  from json import dumps
  open('bravais_summary.json', 'wb').write(dumps(Lfat.as_dict()))
  from cctbx.crystal.crystal_model.serialize import dump_crystal
  for i, subgroup in enumerate(Lfat):
    with open('bravais_setting_%i.json' %(i+1), 'wb') as f:
      dump_crystal(subgroup.refined_crystal, f)
  return


if __name__ == '__main__':
  import sys
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run(sys.argv[1:])
