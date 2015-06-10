from __future__ import division

from libtbx.phil import command_line
from libtbx.utils import Sorry
import iotbx.phil
# from dials.util.command_line import Importer
from dials.util.options import OptionParser
from dials.util.options import flatten_reflections
from dials.util.options import flatten_experiments
from dials.array_family import flex

help_message = '''

This program takes as input the output of dials.index, i.e. experiments.json
and indexed.pickle files. Full refinement of the crystal and experimental
geometry parameters will be performed (by default) in all Bravais settings
that are consistent with the input primitive unit cell. A table is printed
containing various information for each potential Bravais setting, including
the metric fit (a measure of the deviation from the triclinic cell),
the root-mean-square-deviations (rmsd), in mm, between the observed and
predicted spot centroids, the refined unit cell parameters in each Bravais
setting, and the change of basis operator to transform from the triclinic cell
to each Bravais setting.

The program also generates a .json file for each Bravais setting, e.g.
bravais_setting_1.json, which is equivalent to the input experiments.json, but
with the crystal model refined in the chosen Bravais setting. These
bravais_setting_*.json files are suitable as input to dials.refine or
dials.integrate, although the indexed.pickle file will need to be re-indexed
using dials.reindex if the change of basis operator (cb_op) for the chosen
Bravais setting is not the identity operator (a,b,c).

Examples::

  dials.refine_bravais_settings experiments.json indexed.pickle

  dials.refine_bravais_settings experiments.json indexed.pickle nproc=4

'''

phil_scope = iotbx.phil.parse("""
lepage_max_delta = 5
  .type = float
verbosity = 0
  .type = int(value_min=0)
nproc = Auto
  .type = int(value_min=1)
experiment_id = None
  .type = int(value_min=0)
include scope dials.algorithms.refinement.refiner.phil_scope
""", process_includes=True)


def run(args):
  import libtbx.load_env
  usage = "%s experiments.json indexed.pickle [options]" %libtbx.env.dispatcher_name

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_experiments=True,
    read_reflections=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)
  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)
  if len(experiments) == 0:
    parser.print_help()
    return
  elif len(experiments.crystals()) > 1:
    raise Sorry("Only one crystal can be processed at a time: set experiment_id to choose experiment.")

  if params.experiment_id is None:
    params.experiment_id = 0
  assert params.experiment_id < len(experiments)
  experiment = experiments[params.experiment_id]
  assert(len(reflections) == 1)
  reflections = reflections[0]

  from dials.algorithms.indexing.symmetry \
       import refined_settings_factory_from_refined_triclinic

  cb_op_to_primitive = experiment.crystal.get_space_group().info()\
    .change_of_basis_op_to_primitive_setting()
  if experiment.crystal.get_space_group().n_ltr() > 1:
    effective_group = experiment.crystal.get_space_group()\
      .build_derived_reflection_intensity_group(anomalous_flag=True)
    sys_absent_flags = effective_group.is_sys_absent(
      reflections['miller_index'])
    reflections = reflections.select(~sys_absent_flags)
  experiment.crystal.update(experiment.crystal.change_basis(cb_op_to_primitive))
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
    for expt in expts:
      expt.crystal.update(subgroup.refined_crystal)
      expt.detector = subgroup.detector
      expt.beam = subgroup.beam
    dump.experiment_list(
      expts, 'bravais_setting_%i.json' % (int(subgroup.setting_number)))
  return

if __name__ == '__main__':
  import sys
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run(sys.argv[1:])
