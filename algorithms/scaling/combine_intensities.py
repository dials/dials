"""
Optimise the combination of profile and summation intensity values.
"""
from __future__ import print_function
import logging
from dials.array_family import flex
from libtbx.table_utils import simple_table
from cctbx import miller
from dials.algorithms.scaling.scaling_utilities import \
  DialsMergingStatisticsError, calculate_prescaling_correction

logger = logging.getLogger('dials')

logger = logging.getLogger('dials')

def check_for_both_intensities(reflection_tables):
  """Inspect tables to see which have both summation and profile intensities."""
  tables_to_use = []
  tables_to_skip = []
  for i, table in enumerate(reflection_tables):
    assert 'intensity.sum.value' in table
    if 'intensity.prf.value' in table:
      tables_to_use.append(i)
    else:
      tables_to_skip.append(i)
  return tables_to_use, tables_to_skip

def extract_sum_intensity_values(reflection_tables, tables_to_use):
  """Extract the summation intensities from the indicated tables."""
  intensities = flex.double()
  for i in tables_to_use:
    intensities.extend(reflection_tables[i]['intensity.sum.value'].as_double())
  return intensities

def combine_intensities(reflection_table, Imid):
  """Use the given Imid value to perform prf/sum intensity combination."""
  if Imid == 0:
    return _set_intensities_as_prf(reflection_table)
  elif Imid == 1:
    return _set_intensities_as_sum(reflection_table)
  else:
    if 'partiality' in reflection_table:
      i_p = _determine_inverse_partiality(reflection_table)
      sum_intensity = reflection_table['intensity.sum.value'] * i_p
      sum_variance = reflection_table['intensity.sum.variance'] * (i_p ** 2)
    else:
      sum_intensity = reflection_table['intensity.sum.value']
      sum_variance = reflection_table['intensity.sum.variance']
    Icomb, Vcomb = _calculate_combined_raw_intensities(
      reflection_table['intensity.prf.value'], sum_intensity,
      reflection_table['intensity.prf.variance'], sum_variance, Imid)
    conv = calculate_prescaling_correction(reflection_table)
    reflection_table['intensity'] = Icomb * conv
    reflection_table['variance'] = Vcomb * conv * conv
    reflection_table.set_flags(reflection_table['variance'] <= 0.0,
      reflection_table.flags.excluded_for_scaling)
    return reflection_table

def optimise_intensity_combination(reflection_tables, experiment, Imids=None):
  """
  Test combinations of prf/sum intensities to determine optimal Imid value.

  No outlier rejection is performed as it is expected that this function will
  be called after a round of outlier rejection."""

  logger.info("Performing profile/summation intensity optimisation")
  # first analyse which reflection tables have both sum and prf intensities.
  indices_to_use, indices_to_skip = check_for_both_intensities(reflection_tables)
  if len(indices_to_skip) == len(reflection_tables):
    logger.info('''No reflection tables found with both prf and sum values,
no intensity combination can be performed''')
    return 1

  intensities = extract_sum_intensity_values(reflection_tables, indices_to_use)

  # Determine Imid values to try
  if Imids:
    Imid_list = Imids
  else:
    avg = flex.mean(intensities)
    Imid = flex.max(intensities)/10.0
    Imid_list = [0, 1, avg, Imid]
    while (Imid > avg):
      Imid /= 10.0
      Imid_list.append(Imid)

  header = ['Combination', 'CC1/2', 'Rmeas']

  if len(reflection_tables) == 1:
    rows, results = _single_table_combination(
      reflection_tables[0], experiment, Imid_list)
  else:
    rows, results = _multi_table_combination(
      reflection_tables, experiment, Imid_list, indices_to_use)

  st = simple_table(rows, header)
  logger.info(st.format())

  max_key = min(results, key=results.get)
  if max_key == 0:
    logger.info('Profile intensities determined to be best for scaling. \n')
  elif max_key == 1:
    logger.info('Summation intensities determined to be best for scaling. \n')
  else:
    logger.info('Combined intensities with Imid = %s determined to be best for scaling. \n',
      max_key)
  return max_key

def _single_table_combination(reflections, experiment, Imid_list):

  rows = []
  results = {}

  #select good reflections to use
  reflections = _filter_reflections_for_combining(reflections)
  # Calculate prescaling corrections only once
  prescaling_corrections = calculate_prescaling_correction(reflections)
  for Imid in Imid_list:
    Int, Var = _get_Is_from_Imidval(reflections, Imid)
    miller_set = miller.set(crystal_symmetry=experiment.crystal.get_crystal_symmetry(),
      indices=reflections['miller_index'], anomalous_flag=False)
    i_obs = miller.array(miller_set,
      data = Int * prescaling_corrections / reflections['inverse_scale_factor'])
    i_obs.set_observation_type_xray_intensity()
    i_obs.set_sigmas((Var**0.5) * prescaling_corrections / reflections['inverse_scale_factor'])
    array = i_obs.customized_copy(anomalous_flag=False).map_to_asu()

    try:
      merge = array.merge_equivalents()
      rmeas = merge.r_meas()
      cchalf = miller.compute_cc_one_half(unmerged=array)
    except RuntimeError:
      raise DialsMergingStatisticsError("Unable to merge for intensity combination")

    # record the results
    results[Imid] = rmeas
    res_str = {0 : 'prf only', 1 : 'sum only'}
    if not Imid in res_str:
      res_str[Imid] = 'Imid = '+str(round(Imid, 2))
    rows.append([res_str[Imid], str(round(cchalf, 5)), str(round(rmeas, 5))])

  return rows, results

def _multi_table_combination(reflection_tables, experiment, Imid_list, indices_to_use):

  rows = []
  results = {}

  # Calculate prescaling corrections only once
  prescaling_corrections = [None]*len(reflection_tables)
  for i in indices_to_use:
    reflections = reflection_tables[i]
    reflections = _filter_reflections_for_combining(reflections)
    prescaling_corrections[i] = calculate_prescaling_correction(reflections)
  for Imid in Imid_list:
    combined_intensities = flex.double([])
    combined_sigmas = flex.double([])
    combined_scales = flex.double([])
    combined_indices = flex.miller_index([])
    for i in indices_to_use:
      reflections = reflection_tables[i]
      reflections = _filter_reflections_for_combining(reflections)
      Int, Var = _get_Is_from_Imidval(reflections, Imid)
      Int *= prescaling_corrections[i]
      sigma = (Var ** 0.5) * prescaling_corrections[i]
      combined_intensities.extend(Int)
      combined_sigmas.extend(sigma)
      combined_scales.extend(reflections['inverse_scale_factor'])
      combined_indices.extend(reflections['miller_index'])
    # apply scale factor before determining merging stats
    miller_set = miller.set(crystal_symmetry=experiment.crystal.get_crystal_symmetry(),
      indices=combined_indices, anomalous_flag=False)
    i_obs = miller.array(miller_set, data=combined_intensities/combined_scales)
    i_obs.set_observation_type_xray_intensity()
    i_obs.set_sigmas(combined_sigmas/combined_scales)
    array = i_obs.customized_copy(anomalous_flag=False).map_to_asu()
    try:
      merge = array.merge_equivalents()
      rmeas = merge.r_meas()
      cchalf = miller.compute_cc_one_half(unmerged=array)
    except RuntimeError:
      raise DialsMergingStatisticsError("Unable to merge for intensity combination")

    # record the results
    results[Imid] = rmeas
    res_str = {0 : 'prf only', 1 : 'sum only'}
    if not Imid in res_str:
      res_str[Imid] = 'Imid = '+str(round(Imid, 2))
    rows.append([res_str[Imid], str(round(cchalf, 5)), str(round(rmeas, 5))])

  return rows, results



### Helper functions for combine_intensities

def _get_Is_from_Imidval(reflections, Imid):
  """Intepret the Imid value to extract and return the Icomb and Vcomb values."""
  if Imid == 0: #special value to trigger prf
    Int = reflections['intensity.prf.value']
    Var = reflections['intensity.prf.variance']
  elif Imid == 1: #special value to trigger sum
    if 'partiality' in reflections:
      Int = reflections['intensity.sum.value']/reflections['partiality']
      Var = reflections['intensity.sum.variance']/(reflections['partiality']**2)
    else:
      Int = reflections['intensity.sum.value']
      Var = reflections['intensity.sum.variance']
  else:
    if 'partiality' in reflections:
      Int, Var = _calculate_combined_raw_intensities(
        reflections['intensity.prf.value'],
        reflections['intensity.sum.value']/reflections['partiality'],
        reflections['intensity.prf.variance'],
        reflections['intensity.sum.variance']/(reflections['partiality']**2),
        Imid)
    else:
      Int, Var = _calculate_combined_raw_intensities(
        reflections['intensity.prf.value'],
        reflections['intensity.sum.value'],
        reflections['intensity.prf.variance'],
        reflections['intensity.sum.variance'],
        Imid)
  return Int, Var

def _filter_reflections_for_combining(reflections):
  bad_sel = reflections.get_flags(
    reflections.flags.bad_for_scaling, all=False) | \
    (reflections['intensity.prf.variance'] <= 0) | \
    (reflections['intensity.sum.variance'] <= 0) | \
    (reflections['inverse_scale_factor'] <= 0)
  reflections = reflections.select(~bad_sel)
  integrated = reflections.get_flags(
    reflections.flags.integrated, all=True) #might not all have both prf/sum
  reflections = reflections.select(integrated)
  if 'partiality' in reflections:
    reflections = reflections.select(reflections['partiality'] > 0)
  return reflections

def _determine_inverse_partiality(reflections):
  inverse_partiality = flex.double(reflections.size(), 1.0)
  nonzero_partiality_sel = reflections['partiality'] > 0.0
  good_refl = reflections.select(nonzero_partiality_sel)
  inverse_partiality.set_selected(nonzero_partiality_sel.iselection(),
    1.0/good_refl['partiality'])
  return inverse_partiality

def _set_intensities_as_sum(reflections):
  """Set the 'intensity' and 'variance' column to be the summation values, with
  the prescaling corrections applied."""
  conv = calculate_prescaling_correction(reflections)
  reflections['intensity'] = reflections['intensity.sum.value'] * conv
  reflections['variance'] = reflections['intensity.sum.variance'] * conv * conv
  if 'partiality' in reflections:
    inverse_partiality = _determine_inverse_partiality(reflections)
    reflections['intensity'] *= inverse_partiality
    reflections['variance'] *= (inverse_partiality**2)
  reflections.set_flags(reflections['variance'] <= 0.0,
    reflections.flags.excluded_for_scaling)
  return reflections

def _set_intensities_as_prf(reflections):
  conv = calculate_prescaling_correction(reflections)
  reflections['intensity'] = reflections['intensity.prf.value'] * conv
  reflections['variance'] = reflections['intensity.prf.variance'] * conv * conv
  reflections.set_flags(reflections['variance'] <= 0.0,
    reflections.flags.excluded_for_scaling)
  return reflections

def _calculate_combined_raw_intensities(Iprf, Isum, Vprf, Vsum, Imid):
  """Use partiality-corrected Isum, alongside Iprf to calculate
  combined raw intensities."""
  w = 1.0/(1.0 + (Isum/Imid)**3)
  w.set_selected(Isum <= 0, 1.0)
  Icomb = (w * Iprf) + ((1.0 - w) * Isum)
  Vcomb = (w * Vprf) + ((1.0 - w) * Vsum)
  return Icomb, Vcomb


