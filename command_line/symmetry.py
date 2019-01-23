from __future__ import division, absolute_import, print_function

import logging
logger = logging.getLogger('dials.command_line.symmetry')

import copy

from cctbx import crystal
from cctbx import miller
import iotbx.phil

from dials.array_family import flex
from dials.util import log
from dials.util.options import OptionParser
from dials.util.options import flatten_experiments, flatten_reflections
from dials.util.multi_dataset_handling import assign_unique_identifiers,\
  parse_multiple_datasets
from dials.algorithms.symmetry.determine_space_group import determine_space_group
from dials.algorithms.scaling.outlier_rejection import reject_outliers


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
  experiments = "reindexed_experiments.json"
    .type = path
  reflections = "reindexed_reflections.pickle"
    .type = path
  json = dials.symmetry.json
    .type = path
}

''', process_includes=True)

class symmetry(object):
  def __init__(self, experiments, reflections, params=None):
    if params is None:
      params = phil_scope.extract()
    self._params = params

    # transform models into miller arrays
    datasets = self._miller_arrays_from_experiments_reflections(
      experiments, reflections)

    result = determine_space_group(
      datasets, normalisation=self._params.normalisation,
      d_min=self._params.d_min,
      min_i_mean_over_sigma_mean=self._params.min_i_mean_over_sigma_mean,
      relative_length_tolerance=self._params.relative_length_tolerance,
      absolute_angle_tolerance=self._params.absolute_angle_tolerance)
    logger.info(result)

    if params.output.json is not None:
      result.as_json(filename=params.output.json)

    self._export_experiments_reflections(experiments, reflections, result)

  def _miller_arrays_from_experiments_reflections(self,
                                                  experiments, reflections):
    miller_arrays = []

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
          intensity_choice.append('profile')
          intensity_to_use = 'prf'
        else:
          intensity_to_use = 'sum'

      refl = filter_reflection_table(refl, intensity_choice, min_isigi=-5,
        filter_ice_rings=False, combine_partials=True,
        partiality_threshold=self._params.partiality_threshold)
      assert refl.size() > 0
      if intensity_to_use != 'scale':
        try:
          refl['intensity'] = refl['intensity.'+intensity_to_use+'.value']
          refl['variance'] = refl['intensity.'+intensity_to_use+'.variance']
        except RuntimeError:
          intensity_to_use = 'sum'
          refl['intensity'] = refl['intensity.sum.value']
          refl['variance'] = refl['intensity.sum.variance']
        refl = reject_outliers(refl, expt, method='simple', zmax=12.0)
        refl = refl.select(~refl.get_flags(refl.flags.outlier_in_scaling))
      data = refl['intensity.'+intensity_to_use+'.value']
      variances = refl['intensity.'+intensity_to_use+'.variance']

      miller_indices = refl['miller_index']
      assert variances.all_gt(0)
      sigmas = flex.sqrt(variances)

      miller_set = miller.set(crystal_symmetry, miller_indices, anomalous_flag=False)
      intensities = miller.array(miller_set, data=data, sigmas=sigmas)
      intensities.set_observation_type_xray_intensity()
      intensities.set_info(miller.array_info(
        source='DIALS',
        source_type='pickle'
      ))
      miller_arrays.append(intensities)

    return miller_arrays

  def _export_experiments_reflections(self, experiments, reflections, result):
    from dxtbx.serialize import dump
    from rstbx.symmetry.constraints import parameter_reduction
    reindexed_experiments = copy.deepcopy(experiments)
    reindexed_reflections = flex.reflection_table()
    cb_op_inp_best = result.best_solution.subgroup['cb_op_inp_best'] \
      * result.cb_op_inp_min
    best_subsym = result.best_solution.subgroup['best_subsym']
    for i, expt in enumerate(reindexed_experiments):
      expt.crystal = expt.crystal.change_basis(cb_op_inp_best)
      expt.crystal.set_space_group(
        best_subsym.space_group().build_derived_acentric_group())
      S = parameter_reduction.symmetrize_reduce_enlarge(
        expt.crystal.get_space_group())
      S.set_orientation(expt.crystal.get_B())
      S.symmetrize()
      expt.crystal.set_B(S.orientation.reciprocal_matrix())
      reindexed_refl = copy.deepcopy(reflections[i])
      reindexed_refl['miller_index'] = cb_op_inp_best.apply(
        reindexed_refl['miller_index'])
      reindexed_reflections.extend(reindexed_refl)
    logger.info(
      'Saving reindexed experiments to %s' % self._params.output.experiments)
    dump.experiment_list(reindexed_experiments, self._params.output.experiments)
    logger.info('Saving %s reindexed reflections to %s' % (
      len(reindexed_reflections), self._params.output.reflections))
    reindexed_reflections.as_pickle(self._params.output.reflections)


help_message = '''
This program implements the methods of
`POINTLESS <http://www.ccp4.ac.uk/html/pointless.html>`_ (
`Evans, P. (2006). Acta Cryst. D62, 72-82. <https://doi.org/10.1107/S0907444905036693>`_ and
`Evans, P. R. (2011). Acta Cryst. D67, 282-292. <https://doi.org/10.1107/S090744491003982X>`_)
for scoring and determination of Laue group symmetry.

The program takes as input a set of one or more integrated experiments and
reflections.

Examples::

  dials.symmetry experiments.json reflections.pickle

'''


def run(args):
  usage = "dials.cosym [options] experiments.json reflections.pickle"

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_reflections=True,
    read_datablocks=False,
    read_experiments=True,
    check_format=False,
    epilog=help_message
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
  reflections = parse_multiple_datasets(reflections)
  experiments, reflections = assign_unique_identifiers(
    experiments, reflections)

  if len(experiments) == 0 and len(reflections) == 0:
    parser.print_help()
    return

  symmetry(experiments, reflections, params=params)


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
