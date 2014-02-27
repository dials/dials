from __future__ import division

from libtbx.phil import command_line
import iotbx.phil
from dials.util.command_line import Importer

master_phil_scope = iotbx.phil.parse("""
include scope dials.data.refinement.phil_scope
verbosity = 0
  .type = int(value_min=0)
nproc = Auto
  .type = int(value_min=1)
""", process_includes=True)

master_params = master_phil_scope.fetch().extract()


def run(args):
  importer = Importer(args, check_format=False)
  if len(importer.experiments) == 0:
    print "No ExperimentList could be constructed"
    return
  elif len(importer.experiments) > 1:
    raise RuntimeError("Only one Experiment can be processed at a time")
  experiments = importer.experiments
  experiment = experiments[0]
  assert len(importer.reflections) == 1
  reflections = importer.reflections[0]
  args = importer.unhandled_arguments
  cmd_line = command_line.argument_interpreter(master_params=master_phil_scope)
  working_phil = cmd_line.process_and_fetch(args=args)
  working_phil.show()
  params = working_phil.extract()

  from dials.algorithms.indexing.symmetry \
       import refined_settings_factory_from_refined_triclinic

  Lfat = refined_settings_factory_from_refined_triclinic(
    params, experiment, reflections,
    nproc=params.nproc, refiner_verbosity=params.verbosity)
  Lfat.labelit_printout()
  from json import dumps
  open('bravais_summary.json', 'wb').write(dumps(Lfat.as_dict()))
  from dials.model.serialize import dump
  import copy
  for subgroup in Lfat:
    expts = copy.deepcopy(experiments)
    expts[0].crystal = subgroup.refined_crystal
    dump.experiment_list(
      expts, 'bravais_setting_%i.json' % (int(subgroup.setting_number)))
  return

if __name__ == '__main__':
  import sys
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run(sys.argv[1:])
