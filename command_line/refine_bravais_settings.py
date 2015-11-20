from __future__ import division

from cStringIO import StringIO
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
crystal_id = None
  .type = int(value_min=0)
normalise = False
  .type = bool
  .help = "Normalise intensities before calculating correlation coefficients."
normalise_bins = 0
  .type = int
  .help = "Number of resolution bins for normalisation"

output {
  directory = None
    .type = str
  log = dials.refine_bravais_settings.log
    .type = str
  debug_log = dials.refine_bravais_settings.debug.log
    .type = str
}

include scope dials.algorithms.refinement.refiner.phil_scope
""", process_includes=True)


def run(args):
  from dials.util import log
  from logging import info
  import libtbx.load_env
  usage = "%s experiments.json indexed.pickle [options]" %libtbx.env.dispatcher_name

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_experiments=True,
    read_reflections=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=False)

  # Configure the logging
  log.config(info=params.output.log, debug=params.output.debug.log)

  from dials.util.version import dials_version
  info(dials_version())

  # Log the diff phil
  diff_phil = parser.diff_phil.as_str()
  if diff_phil is not '':
    info('The following parameters have been modified:\n')
    info(diff_phil)

  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)
  assert(len(reflections) == 1)
  reflections = reflections[0]

  if len(experiments) == 0:
    parser.print_help()
    return
  elif len(experiments.crystals()) > 1:
    if params.crystal_id is not None:
      assert params.crystal_id < len(experiments.crystals())
      experiment_ids = experiments.where(crystal=experiments.crystals()[params.crystal_id])
      from dxtbx.model.experiment.experiment_list import ExperimentList
      experiments = ExperimentList([experiments[i] for i in experiment_ids])
      refl_selections = [reflections['id'] == i for i in experiment_ids]
      for i, sel in enumerate(refl_selections):
        reflections['id'].set_selected(sel, i)
    else:
      raise Sorry("Only one crystal can be processed at a time: set crystal_id to choose experiment.")

  from dials.algorithms.indexing.symmetry \
       import refined_settings_factory_from_refined_triclinic

  cb_op_to_primitive = experiments[0].crystal.get_space_group().info()\
    .change_of_basis_op_to_primitive_setting()
  if experiments[0].crystal.get_space_group().n_ltr() > 1:
    effective_group = experiments[0].crystal.get_space_group()\
      .build_derived_reflection_intensity_group(anomalous_flag=True)
    sys_absent_flags = effective_group.is_sys_absent(
      reflections['miller_index'])
    reflections = reflections.select(~sys_absent_flags)
  experiments[0].crystal.update(experiments[0].crystal.change_basis(cb_op_to_primitive))
  miller_indices = reflections['miller_index']
  miller_indices = cb_op_to_primitive.apply(miller_indices)
  reflections['miller_index'] = miller_indices

  Lfat = refined_settings_factory_from_refined_triclinic(
    params, experiments, reflections, lepage_max_delta=params.lepage_max_delta,
    nproc=params.nproc, refiner_verbosity=params.verbosity)
  s = StringIO()
  Lfat.labelit_printout(out=s)
  info(s.getvalue())
  from json import dumps
  from os.path import join
  open(join(params.output.directory, 'bravais_summary.json'), 'wb').write(dumps(Lfat.as_dict()))
  from dxtbx.serialize import dump
  import copy
  for subgroup in Lfat:
    expts = copy.deepcopy(experiments)
    for expt in expts:
      expt.crystal.update(subgroup.refined_crystal)
      expt.detector = subgroup.detector
      expt.beam = subgroup.beam
    dump.experiment_list(
      expts, join(params.output.directory, 'bravais_setting_%i.json' % (int(subgroup.setting_number))))
  return

if __name__ == '__main__':
  import sys
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run(sys.argv[1:])
