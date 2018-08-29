from __future__ import division, absolute_import, print_function

import logging
logger = logging.getLogger('dials.command_line.symmetry')

import copy
import os

from cctbx import crystal
from cctbx import miller
import iotbx.phil
from iotbx.reflection_file_reader import any_reflection_file

from dials.array_family import flex
from dials.util import log
from dials.util.options import OptionParser
from dials.util.options import flatten_experiments, flatten_reflections
from dials.algorithms.symmetry.determine_space_group import determine_space_group


phil_scope = iotbx.phil.parse('''\
d_min = Auto
  .type = float(value_min=0)

min_i_mean_over_sigma_mean = 4
  .type = float(value_min=0)

min_cc_half = 0.6
  .type = float(value_min=0, value_max=1)

batch = None
  .type = ints(value_min=0, size=2)

normalisation = kernel quasi ml_iso *ml_aniso
  .type = choice

lattice_group = None
  .type = space_group

verbosity = 1
  .type = int(value_min=0)
  .help = "The verbosity level"

seed = 230
  .type = int(value_min=0)

relative_length_tolerance = 0.05
  .type = float(value_min=0)

absolute_angle_tolerance = 2
  .type = float(value_min=0)

partiality_threshold = 0.99
  .type = float
  .help = "Use only reflections with a partiality above this threshold."

output {
  log = dials.symmetry.log
    .type = str
  debug_log = dials.symmetry.debug.log
    .type = str
  suffix = "_reindexed"
    .type = str
  experiments = "reindexed_experiments.json"
    .type = path
  reflections = "reindexed_reflections.pickle"
    .type = path
  json = dials.symmetry.json
    .type = path
}

''', process_includes=True)

def run(args):

  parser = OptionParser(
    #usage=usage,
    phil=phil_scope,
    read_reflections=True,
    read_datablocks=False,
    read_experiments=True,
    check_format=False,
    #epilog=help_message
  )

  params, options, args = parser.parse_args(
    args=args, show_diff_phil=False, return_unhandled=True)

  # Configure the logging
  log.config(
    params.verbosity,
    info=params.output.log,
    debug=params.output.debug_log)

  from dials.util.version import dials_version
  logger.info(dials_version())

  # Log the diff phil
  diff_phil = parser.diff_phil.as_str()
  if diff_phil is not '':
    logger.info('The following parameters have been modified:\n')
    logger.info(diff_phil)

  if params.seed is not None:
    import random
    flex.set_random_seed(params.seed)
    random.seed(params.seed)

  datasets = []

  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)

  if len(experiments) or len(reflections):
    if len(reflections) == 1:
      reflections_input = reflections[0]
      reflections = []
      for i in range(len(experiments)):
        reflections.append(reflections_input.select(reflections_input['id'] == i))

    if len(experiments) > len(reflections):
      flattened_reflections = []
      for refl in reflections:
        for i in range(0, flex.max(refl['id'])+1):
          sel = refl['id'] == i
          flattened_reflections.append(refl.select(sel))
      reflections = flattened_reflections

    assert len(experiments) == len(reflections)

    for expt, refl in zip(experiments, reflections):
      crystal_symmetry = crystal.symmetry(
        unit_cell=expt.crystal.get_unit_cell(),
        space_group=expt.crystal.get_space_group())

      from dials.util.filter_reflections import filter_reflection_table
      if 'intensity.scale.value' in refl:
        intensity_choice = ['scale']
        intensity_to_use = 'scale'
      else:
        assert 'intensity.sum.value' in refl
        intensity_choice = ['sum']
        if 'intensity.prf.value' in refl:
          intensity_choice.append('prf')
          intensity_to_use = 'prf'
        else:
          intensity_to_use = 'sum'

      refl = filter_reflection_table(refl, intensity_choice, min_isigi=-5,
        filter_ice_rings=False, combine_partials=True,
        partiality_threshold=params.partiality_threshold)
      data = refl['intensity.'+intensity_to_use+'.value']
      variances = refl['intensity.'+intensity_to_use+'.variance']

      miller_indices = refl['miller_index']
      assert variances.all_gt(0)
      sigmas = flex.sqrt(variances)

      miller_set = miller.set(
        crystal_symmetry, miller_indices, anomalous_flag=True)
      intensities = miller.array(miller_set, data=data, sigmas=sigmas)
      intensities.set_observation_type_xray_intensity()
      datasets.append(intensities)

  files = args
  for file_name in files:
    reader = any_reflection_file(file_name)
    assert reader.file_type() == 'ccp4_mtz'

    as_miller_arrays = reader.as_miller_arrays(merge_equivalents=False)
    intensities_prf = [ma for ma in as_miller_arrays
                       if ma.info().labels == ['IPR', 'SIGIPR']]
    intensities_sum = [ma for ma in as_miller_arrays
                       if ma.info().labels == ['I', 'SIGI']]
    if len(intensities_prf):
      intensities = intensities_prf[0]
    else:
      assert len(intensities_sum), 'No intensities found in input file.'
      intensities = intensities_sum[0]
    batches = [ma for ma in as_miller_arrays
               if ma.info().labels == ['BATCH']]
    if len(batches):
      batches = batches[0]
    else:
      batches = None
    mtz_object = reader.file_content()
    intensities = intensities.customized_copy(
      anomalous_flag=True,
      indices=mtz_object.extract_original_index_miller_indices())

    intensities.set_observation_type_xray_intensity()
    if params.batch is not None:
      assert batches is not None
      bmin, bmax = params.batch
      assert bmax >= bmin
      sel = (batches.data() >= bmin) & (batches.data() <= bmax)
      assert sel.count(True) > 0
      intensities = intensities.select(sel)

    datasets.append(intensities)

  result = determine_space_group(
    datasets, normalisation=params.normalisation,
    d_min=params.d_min,
    min_i_mean_over_sigma_mean=params.min_i_mean_over_sigma_mean,
    relative_length_tolerance=params.relative_length_tolerance,
    absolute_angle_tolerance=params.absolute_angle_tolerance)

  if params.output.json is not None:
    result.as_json(filename=params.output.json)

  if (len(experiments) and len(reflections) and
      params.output.reflections is not None and
      params.output.experiments is not None):
    from dxtbx.serialize import dump
    from rstbx.symmetry.constraints import parameter_reduction
    reindexed_experiments = copy.deepcopy(experiments)
    reindexed_reflections = flex.reflection_table()
    cb_op_inp_best = result.best_solution.subgroup['cb_op_inp_best'] * result.cb_op_inp_min
    best_subsym = result.best_solution.subgroup['best_subsym']
    for i, expt in enumerate(reindexed_experiments):
      expt.crystal = expt.crystal.change_basis(cb_op_inp_best)
      expt.crystal.set_space_group(best_subsym.space_group().build_derived_acentric_group())
      S = parameter_reduction.symmetrize_reduce_enlarge(expt.crystal.get_space_group())
      S.set_orientation(expt.crystal.get_B())
      S.symmetrize()
      expt.crystal.set_B(S.orientation.reciprocal_matrix())
      reindexed_refl = copy.deepcopy(reflections[i])
      reindexed_refl['miller_index'] = cb_op_inp_best.apply(
        reindexed_refl['miller_index'])
      reindexed_refl['id'] = flex.int(len(reindexed_refl), i)
      reindexed_reflections.extend(reindexed_refl)
    logger.info('Saving reindexed experiments to %s' % params.output.experiments)
    dump.experiment_list(reindexed_experiments, params.output.experiments)
    logger.info('Saving reindexed reflections to %s' % params.output.reflections)
    reindexed_reflections.as_pickle(params.output.reflections)

  elif params.output.suffix is not None:
    cb_op_inp_best = result.best_solution.subgroup['cb_op_inp_best'] * result.cb_op_inp_min
    best_subsym = result.best_solution.subgroup['best_subsym']
    space_group = best_subsym.space_group().build_derived_acentric_group()
    for file_name in files:
      basename = os.path.basename(file_name)
      out_name = os.path.splitext(basename)[0] + params.output.suffix + ".mtz"
      reader = any_reflection_file(file_name)
      assert reader.file_type() == 'ccp4_mtz'
      mtz_object = reader.file_content()
      if not cb_op_inp_best.is_identity_op():
        mtz_object.change_basis_in_place(cb_op_inp_best)
      mtz_object.set_space_group_info(space_group.info())
      for cryst in mtz_object.crystals():
        cryst.set_unit_cell_parameters(best_subsym.unit_cell().parameters())
      mtz_object.write(out_name)
      logger.info('Saving reindexed reflections to %s' % out_name)


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
