from __future__ import division

from libtbx.phil import command_line
from libtbx.utils import Sorry
import iotbx.phil
from dials.util.command_line import Importer
from dials.array_family import flex

master_phil_scope = iotbx.phil.parse("""
include scope dials.data.refinement.phil_scope
lepage_max_delta = 5
  .type = float
verbosity = 0
  .type = int(value_min=0)
nproc = Auto
  .type = int(value_min=1)
experiment_id = None
  .type = int(value_min=0)
""", process_includes=True)

master_params = master_phil_scope.fetch().extract()


def run(args):
  if len(args) == 0:
    from libtbx.utils import Usage
    import libtbx.load_env
    from cStringIO import StringIO
    usage_message = """\
%s experiments.json indexed.pickle [options]

Parameters:
""" %libtbx.env.dispatcher_name
    s = StringIO()
    master_phil_scope.show(out=s)
    usage_message += s.getvalue()
    raise Usage(usage_message)
  importer = Importer(args, check_format=False)
  args = importer.unhandled_arguments
  cmd_line = command_line.argument_interpreter(master_params=master_phil_scope)
  working_phil = cmd_line.process_and_fetch(args=args)
  working_phil.show()
  params = working_phil.extract()
  if len(importer.experiments) == 0:
    print "No ExperimentList could be constructed"
    return
  elif len(importer.experiments) > 1 and params.experiment_id is None:
    raise Sorry("More than one experiment present: set experiment_id to choose experiment.")
  experiments = importer.experiments
  if params.experiment_id is None:
    params.experiment_id = 0
  assert params.experiment_id < len(experiments)
  experiment = experiments[params.experiment_id]
  assert len(importer.reflections) == 1
  reflections = importer.reflections[0]

  from dials.algorithms.indexing.symmetry \
       import refined_settings_factory_from_refined_triclinic

  from dxtbx.model.crystal import crystal_model
  cb_op_to_primitive = experiment.crystal.get_space_group().info()\
    .change_of_basis_op_to_primitive_setting()
  if experiment.crystal.get_space_group().n_ltr() > 1:
    effective_group = experiment.crystal.get_space_group()\
      .build_derived_reflection_intensity_group(anomalous_flag=True)
    sys_absent_flags = effective_group.is_sys_absent(
      reflections['miller_index'])
    reflections = reflections.select(~sys_absent_flags)
  experiment.crystal = experiment.crystal.change_basis(cb_op_to_primitive)
  reflections = reflections.select(reflections['id'] == params.experiment_id)
  miller_indices = reflections['miller_index']
  miller_indices = cb_op_to_primitive.apply(miller_indices)
  reflections['miller_index'] = miller_indices
  reflections['id'] = flex.int(len(reflections), 0)

  Lfat = refined_settings_factory_from_refined_triclinic(
    params, experiment, reflections, lepage_max_delta=params.lepage_max_delta,
    nproc=params.nproc, refiner_verbosity=params.verbosity)
  Lfat.labelit_printout()
  from json import dumps
  open('bravais_summary.json', 'wb').write(dumps(Lfat.as_dict()))
  from dxtbx.serialize import dump
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
