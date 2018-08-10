"""
Module of library functions, to perform core scaling operations on reflection
tables and experiment lists. Some functions, such as create_scaling_model and
merging statistics calculations are called from the main dials.scale script,
whereas others are provided as library functions for calling from custom
scripts. The functions defined here should ideally only require reflection
tables and ExperimentList objects (and sometimes phil_scope objects if
necessary), and return common dials objects such as reflection tables and
ExperimentLists.
"""
from copy import deepcopy
import logging
import pkg_resources
from libtbx import phil
from mock import Mock
import iotbx.merging_statistics
from iotbx import cif, mtz
from cctbx import miller, crystal
from dials.array_family import flex
from dials.util.options import OptionParser
from dials.algorithms.scaling.model.scaling_model_factory import \
  KBSMFactory
from dials.algorithms.scaling.Ih_table import IhTable
from dials.algorithms.scaling.scaling_utilities import get_next_unique_id, \
  calculate_prescaling_correction

logger = logging.getLogger('dials')

def choose_scaling_intensities(reflection_table, intensity_choice='prf'):
  """Choose which intensities to use for scaling. The LP, QE and
  partiality corrections are also applied. Two new columns are
  added to the reflection table 'intensity' and 'variance', which have
  all corrections applied except an inverse scale factor."""
  conv = calculate_prescaling_correction(reflection_table)
  intstr = 'intensity.'+intensity_choice+'.value'
  if not intstr in reflection_table:
  #Can't find selection, try to choose prf, if not then sum (also catches combine
  # which should not be used at this point)
    if 'intensity.prf.value' in reflection_table:
      intstr = 'intensity.prf.value'
    else:
      assert 'intensity.sum.value' in reflection_table, '''No recognised
        intensity values found.'''
      intstr = 'intensity.sum.value'
    varstr = intstr.rstrip('value') + 'variance'
  else:
    varstr = intstr.rstrip('value') + 'variance'
    logger.info(('{0} intensities will be used for scaling (and mtz \n'
        'output if applicable). \n').format(intstr))

  #prf partial intensities are the 'full' intensity values but sum are not
  if 'partiality' in reflection_table and intstr == 'intensity.sum.value':
    inverse_partiality = flex.double(reflection_table.size(), 1.0)
    nonzero_partiality_sel = reflection_table['partiality'] > 0.0
    good_refl = reflection_table.select(reflection_table['partiality'] > 0.0)
    inverse_partiality.set_selected(nonzero_partiality_sel.iselection(),
      1.0/good_refl['partiality'])
    conv *= inverse_partiality

  reflection_table['intensity'] = reflection_table[intstr] * conv
  reflection_table['variance'] = reflection_table[varstr] * conv * conv

  variance_mask = reflection_table['variance'] <= 0.0
  reflection_table.set_flags(variance_mask,
    reflection_table.flags.excluded_for_scaling)
  return reflection_table

def scale_against_target(reflection_table, experiment, target_reflection_table,
  target_experiment, params=None, model='KB'):
  """Determine scale factors for a single dataset, by scaling against a target
  reflection table. Requires a single reflection table for the reflections to
  scale and the target dataset, and an ExperimentList for both datasets. The
  params option can also be specified, if None then the default scaling
  configuration is used. The scaling model can be specified individually.

  Returns the reflection table, with added columns 'inverse_scale_factor' and
  'inverse_scale_factor_variance'."""

  assert model in ['physical', 'array', 'KB']
  if not params:
    phil_scope = phil.parse('''
      include scope dials.algorithms.scaling.scaling_options.phil_scope
      include scope dials.algorithms.scaling.scaling_refiner.scaling_refinery_phil_scope
    ''', process_includes=True)
    optionparser = OptionParser(phil=phil_scope, check_format=False)
    params, _ = optionparser.parse_args(args=[], quick_parse=True)
    params.__inject__('model', model)

  from dials.algorithms.scaling.scaler_factory import TargetScalerFactory

  reflections = [reflection_table, target_reflection_table]
  experiment.append(target_experiment[0])
  experiments = create_scaling_model(params, experiment, reflections)
  scaler = TargetScalerFactory.create(params, experiments, reflections,
    is_scaled_list=[False, True])
  scaler.perform_scaling()
  scaler.expand_scales_to_all_reflections(calc_cov=True)
  return scaler.unscaled_scalers[0].reflection_table

def scale_single_dataset(reflection_table, experiment, params=None,
    model='physical'):
  """Determine scale factors for a single dataset. Requires a reflection table
  and an ExperimentList with a single experiment. A custom params option can be
  specified, if not the default scaling params option will be used, with default
  configuration options. The model can be individually specified.

  Returns the reflection table, with added columns 'inverse_scale_factor' and
  'inverse_scale_factor_variance'."""

  assert model in ['physical', 'array']
  if not params:
    phil_scope = phil.parse('''
      include scope dials.algorithms.scaling.scaling_options.phil_scope
      include scope dials.algorithms.scaling.scaling_refiner.scaling_refinery_phil_scope
    ''', process_includes=True)
    optionparser = OptionParser(phil=phil_scope, check_format=False)
    params, _ = optionparser.parse_args(args=[], quick_parse=True)
  params.__inject__('model', model)

  from dials.algorithms.scaling.scaler_factory import SingleScalerFactory

  experiments = create_scaling_model(params, experiment, [reflection_table])
  scaler = SingleScalerFactory.create(params, experiments[0], reflection_table)
  scaler.perform_scaling()
  scaler.outlier_rejection_routine()
  scaler.perform_scaling(engine=params.scaling_refinery.full_matrix_engine,
    max_iterations=params.scaling_refinery.full_matrix_max_iterations)
  scaler.expand_scales_to_all_reflections(calc_cov=True)
  scaler.round_of_outlier_rejection()
  return scaler.reflection_table

def create_scaling_model(params, experiments, reflections):
  """Create or load a scaling model for multiple datasets."""
  for i, (exp, refl) in enumerate(zip(experiments, reflections)):
    model = experiments.scaling_models()[i]
    '''if params.scaling_options.target_intensities and i == len(reflections)-1:
      for entry_point in pkg_resources.iter_entry_points('dxtbx.scaling_model_ext'):
        if entry_point.name == 'KB':
          #finds relevant extension in dials.extensions.scaling_model_ext
          factory = entry_point.load().factory()
          exp.scaling_model = factory.create(params, exp, refl)
          exp.scaling_model.set_scaling_model_as_scaled()'''
    if model is not None:
      exp.scaling_model = model
    else:
      for entry_point in pkg_resources.iter_entry_points('dxtbx.scaling_model_ext'):
        if entry_point.name == params.model:
          #finds relevant extension in dials.extensions.scaling_model_ext
          factory = entry_point.load().factory()
          exp.scaling_model = factory.create(params, exp, refl)
  return experiments

def create_Ih_table(experiments, reflections, selections=None, n_blocks=1,
  weighting_scheme=None):
  """Create an Ih table from a list of experiments and reflections. Optionally,
  a selection list can also be given, to select data from each reflection table.
  Allow an unequal number of experiments and reflections, as only need to
  extract one space group value (can optionally check all same if many)."""
  if selections:
    assert len(selections) == len(reflections), """Must have an equal number of
    reflection tables and selections in the input lists."""
  space_group_0 = experiments[0].crystal.get_space_group()
  for experiment in experiments:
    assert experiment.crystal.get_space_group() == space_group_0, """The space
    groups of all experiments must be equal."""
  refl_and_sel_list = []
  for i, reflection in enumerate(reflections):
    if not 'inverse_scale_factor' in reflection:
      reflection['inverse_scale_factor'] = flex.double(reflection.size(), 1.0)
    if selections:
      refl_and_sel_list.append((reflection, selections[i]))
    else:
      refl_and_sel_list.append((reflection, None))
  Ih_table = IhTable(refl_and_sel_list, space_group_0, n_blocks, weighting_scheme)
  return Ih_table

def calculate_merging_statistics(reflection_table, experiments, use_internal_variance):
  """Calculate merging statistics for scaled datasets. Datasets are selected
  from the reflection table based on their id, and a list of dataset statistics
  objects and dataset ids are returned."""
  results = []
  ids = []
  dataset_ids = list(set(reflection_table['id']))
  if len(dataset_ids) == 1:
    results.append(calculate_single_merging_stats(reflection_table,
      experiments[0], use_internal_variance))
    ids.append(dataset_ids[0])
  else:
    for dataset_id in dataset_ids:
      refls = reflection_table.select(reflection_table['id'] == dataset_id)
      results.append(calculate_single_merging_stats(refls, experiments[0],
        use_internal_variance))
      ids.append(dataset_id)
  return results, ids

def calculate_single_merging_stats(reflection_table, experiment,
  use_internal_variance, n_bins=20):
  """Calculate the merging stats for a single dataset."""
  bad_refl_sel = reflection_table.get_flags(
    reflection_table.flags.bad_for_scaling, all=False)
  r_t = reflection_table.select(~bad_refl_sel)
  miller_set = miller.set(
    crystal_symmetry=experiment.crystal.get_crystal_symmetry(),
    indices=r_t['miller_index'], anomalous_flag=False)
  i_obs = miller.array(
    miller_set, data=r_t['intensity']/r_t['inverse_scale_factor'])
  i_obs.set_observation_type_xray_intensity()
  i_obs.set_sigmas((r_t['variance']**0.5)/r_t['inverse_scale_factor'])
  i_obs.set_info(
    miller.array_info(source='DIALS', source_type='reflection_tables'))
  result = iotbx.merging_statistics.dataset_statistics(
    i_obs=i_obs, n_bins=n_bins, anomalous=False, sigma_filtering=None,
    use_internal_variance=use_internal_variance, eliminate_sys_absent=False)
  return result

def intensity_array_from_cif_file(cif_file):
  """Return an intensity miller array from a cif file."""
  model = cif.reader(file_path=cif_file).build_crystal_structures()['1']
  ic = model.structure_factors(anomalous_flag=True, d_min=0.4,
    algorithm='direct').f_calc().as_intensity_array()
  return ic

def create_datastructures_for_target_mtz(experiments, mtz_file):
  """Read a merged mtz file and extract miller indices, intensities and
  variances."""
  m = mtz.object(mtz_file)
  ind = m.extract_miller_indices()
  cols = m.columns()
  col_dict = {c.label() : c for c in cols}
  r_t = flex.reflection_table()
  if 'I' in col_dict: #nice and simple
    r_t['miller_index'] = ind
    r_t['intensity'] = col_dict['I'].extract_values().as_double()
    r_t['variance'] = col_dict['SIGI'].extract_values().as_double()
  elif 'IMEAN' in col_dict: #nice and simple
    r_t['miller_index'] = ind
    r_t['intensity'] = col_dict['IMEAN'].extract_values().as_double()
    r_t['variance'] = col_dict['SIGIMEAN'].extract_values().as_double()
  elif 'I(+)' in col_dict: #need to combine I+ and I- together into target Ih
    if col_dict['I(+)'].n_valid_values() == 0:#use I(-)
      r_t['miller_index'] = ind
      r_t['intensity'] = col_dict['I(-)'].extract_values().as_double()
      r_t['variance'] = col_dict['SIGI(-)'].extract_values().as_double()
    elif col_dict['I(-)'].n_valid_values() == 0:#use I(+)
      r_t['miller_index'] = ind
      r_t['intensity'] = col_dict['I(+)'].extract_values().as_double()
      r_t['variance'] = col_dict['SIGI(+)'].extract_values().as_double()
    else: #Combine both - add together then use Ih table to calculate I and sigma
      r_tplus = flex.reflection_table()
      r_tminus = flex.reflection_table()
      r_tplus['miller_index'] = ind
      r_tplus['intensity'] = col_dict['I(+)'].extract_values().as_double()
      r_tplus['variance'] = col_dict['SIGI(+)'].extract_values().as_double()
      r_tminus['miller_index'] = ind
      r_tminus['intensity'] = col_dict['I(-)'].extract_values().as_double()
      r_tminus['variance'] = col_dict['SIGI(-)'].extract_values().as_double()
      r_tplus.extend(r_tminus)
      r_tplus.set_flags(flex.bool(r_tplus.size(), False),
        r_tplus.flags.bad_for_scaling)
      r_tplus = r_tplus.select(r_tplus['variance'] != 0.0)
      Ih_table = create_Ih_table([experiments[0]], [r_tplus]).blocked_data_list[0]
      r_t['intensity'] = Ih_table.Ih_values
      inv_var = (Ih_table.weights * Ih_table.h_index_matrix) * Ih_table.h_expand_matrix
      r_t['variance'] = 1.0/inv_var
      r_t['miller_index'] = Ih_table.miller_index
  else:
    assert 0, """Unrecognised intensities in mtz file."""

  r_t.set_flags(flex.bool(r_t.size(), True), r_t.flags.integrated)
  exp = deepcopy(experiments[0]) #copy exp for space group -
    #any other necessary reason or can this attribute be added?
  used_ids = experiments.identifiers()
  unique_id = get_next_unique_id(0, used_ids)
  exp.identifier = str(unique_id)
  r_t.experiment_identifiers()[unique_id] = str(unique_id)

  # create a new KB scaling model for the target and set as scaled to fix scale
  # for targeted scaling.
  params = Mock()
  params.parameterisation.decay_term.return_value = False
  params.parameterisation.scale_term.return_value = True
  exp.scaling_model = KBSMFactory.create(params, [], [])
  exp.scaling_model.set_scaling_model_as_scaled() #Set as scaled to fix scale.

  return exp, r_t

def create_datastructures_for_structural_model(reflections, experiments,
    cif_file):
  """Read a cif file, calculate intensities and scale them to the average
  intensity of the reflections. Return an experiment and reflection table to
  be used for the structural model in scaling."""

  # read model, compute Fc, square to F^2
  ic = intensity_array_from_cif_file(cif_file)
  exp = deepcopy(experiments[0])
  params = Mock()
  params.parameterisation.decay_term.return_value = False
  params.parameterisation.scale_term.return_value = True
  exp.scaling_model = KBSMFactory.create(params, [], [])
  exp.scaling_model.set_scaling_model_as_scaled() #Set as scaled to fix scale.

  # Now put the calculated I's on roughly a common scale with the data.
  miller_indices = flex.miller_index([])
  intensities = flex.double([])

  for refl in reflections:
    miller_indices.extend(refl['miller_index'])
    intensities.extend(refl['intensity.prf.value'])
  miller_set = miller.set(crystal_symmetry=crystal.symmetry(
    space_group=experiments[0].crystal.get_space_group()),
    indices=miller_indices, anomalous_flag=True)
  idata = miller.array(miller_set, data=intensities)

  match = idata.match_indices(ic)
  pairs = match.pairs()

  icalc = flex.double()
  iobs = flex.double()
  miller_idx = flex.miller_index()
  for p in pairs:
    # Note : will create miller_idx duplicates in i_calc - problem?
    iobs.append(idata.data()[p[0]])
    icalc.append(ic.data()[p[1]])
    miller_idx.append(ic.indices()[p[1]])

  icalc *= flex.sum(iobs) / flex.sum(icalc)

  rt = flex.reflection_table()
  rt['intensity'] = icalc
  rt['miller_index'] = miller_idx

  return exp, rt
