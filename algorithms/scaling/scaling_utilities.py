"""
Module of utility functions for scaling.
"""

from __future__ import print_function
import sys
import logging
from time import time
import copy
from math import pi, acos
from dials.array_family import flex
from libtbx.table_utils import simple_table
from libtbx.utils import Sorry
import iotbx.merging_statistics
from cctbx import miller
from cctbx import uctbx
from dxtbx.model.experiment_list import ExperimentListDumper, ExperimentList
from dials_scaling_ext import create_sph_harm_table, calc_theta_phi,\
  rotate_vectors_about_axis


logger = logging.getLogger('dials')

try:
  import resource
  import platform
  def log_memory_usage():
    # getrusage returns kb on linux, bytes on mac
    units_per_mb = 1024
    if platform.system() == "Darwin":
      units_per_mb = 1024*1024
    logger.debug('Memory usage: %.1f MB' % (int(resource.getrusage(
      resource.RUSAGE_SELF).ru_maxrss) / units_per_mb))
except ImportError:
  def log_memory_usage():
    pass

class DialsMergingStatisticsError(Exception):
  """Raised when iotbx merging statistics fails."""
  pass

class BadDatasetForScalingException(Exception):
  """Raised when a selection leaves no further good reflections."""
  pass

class Reasons(object):

  def __init__(self):
    self.reasons = {}

  def add_reason(self, text, number):
    self.reasons[text] = number

  def __repr__(self):
    reasonlist = ['criterion: %s, reflections: %s\n' % (k, v) for (k, v) in self.reasons.iteritems() if v > 0]
    return 'Reflections passing individual criteria:\n'+''.join(reasonlist)

def save_experiments(experiments, filename):
  """Save the experiments json."""
  st = time()
  logger.info('Saving the experiments to %s', filename)
  dump = ExperimentListDumper(experiments)
  with open(filename, "w") as outfile:
    outfile.write(dump.as_json(split=True))
  logger.info('Time taken: %g', (time() - st))

def save_reflections(reflection_table, filename):
  """Save the scaled reflections."""
  st = time()
  logger.info('Saving the scaled reflections to %s', filename)
  reflection_table.as_pickle(filename)
  logger.info('Time taken: %g', (time() - st))

'''def calc_sigmasq(jacobian_transpose, var_cov):
  sigmasq = flex.float([])
  #note: must be a faster way to do this next bit? - in c++?
  for col in jacobian_transpose.cols(): #iterating over reflections
    a = flex.double(col.as_dense_vector())
    var = (a * var_cov) * a
    sigmasq.append(flex.sum(var))
  return sigmasq.as_double()'''

def calc_crystal_frame_vectors(reflection_table, experiments):
  """Calculate the diffraction vectors in the crystal frame."""
  reflection_table['s0'] = flex.vec3_double(
    [experiments.beam.get_s0()]*len(reflection_table))
  rot_axis = flex.vec3_double([experiments.goniometer.get_rotation_axis()])
  angles = reflection_table['phi'] * -1.0 * pi / 180 #want to do an inverse rot.
  reflection_table['s1c'] = rotate_vectors_about_axis(rot_axis,
    reflection_table['s1'], angles)
  reflection_table['s0c'] = rotate_vectors_about_axis(rot_axis,
    reflection_table['s0'], angles)
  reflection_table['s1c'] = align_rotation_axis_along_z(rot_axis,
    reflection_table['s1c'])
  reflection_table['s0c'] = align_rotation_axis_along_z(rot_axis,
    reflection_table['s0c'])
  return reflection_table

def align_rotation_axis_along_z(exp_rot_axis, vectors):
  """Rotate the coordinate system such that the exp_rot_axis is along z."""
  if list(exp_rot_axis) == [(0.0, 0.0, 1.0)]:
    return vectors
  (ux, uy, uz) = exp_rot_axis[0][0], exp_rot_axis[0][1], exp_rot_axis[0][2]
  cross_prod_uz = flex.vec3_double([(uy, -1.0*ux, 0.0)])
  angle_between_u_z = +1.0 * acos(uz/((ux**2 + uy**2 + uz**2)**0.5))
  phi = flex.double(vectors.size(), angle_between_u_z)
  new_vectors = rotate_vectors_about_axis(cross_prod_uz, vectors, phi)
  return flex.vec3_double(new_vectors)

def sph_harm_table(reflection_table, experiments, lmax):
  """Calculate the spherical harmonic table for a spherical
    harmonic absorption correction."""
  reflection_table = calc_crystal_frame_vectors(reflection_table, experiments)
  theta_phi = calc_theta_phi(reflection_table['s0c'])
  theta_phi_2 = calc_theta_phi(reflection_table['s1c'])
  sph_h_t = create_sph_harm_table(theta_phi, theta_phi_2, lmax)
  return sph_h_t

def quasi_normalisation(reflection_table, experiment):
  """Calculate normalised intensity (Esq) values for reflections, for the purpose
  of selecting subsets based on Esq for scaling. If more involved analyses of
  normalised intensities are needed, then it may be necessary to split this
  procedure to handle acentric and centric reflections separately."""

  logger.info('Calculating normalised intensity values to select a reflection \n'
    'subset for scaling. \n')
  logger.debug('Negative intensities are set to zero for the purpose of \n'
    'calculating mean intensity values for resolution bins. This is to avoid \n'
    'spuriously high E^2 values due to a mean close to zero and should only \n'
    'affect the E^2 values of the highest resolution bins. \n')

  good_refl_sel = ~reflection_table.get_flags(
    reflection_table.flags.bad_for_scaling, all=False)
  rt_subset = reflection_table.select(good_refl_sel)

  # Scaling subset is data that has not been flagged as bad or excluded
  miller_set = miller.set(
    crystal_symmetry=experiment.crystal.get_crystal_symmetry(),
    indices=rt_subset['miller_index'])

  #handle negative reflections to minimise effect on mean I values.
  rt_subset['intensity_for_norm'] = copy.deepcopy(rt_subset['intensity'])
  rt_subset['intensity_for_norm'].set_selected(rt_subset['intensity'] < 0.0, 0.0)
  miller_array = miller.array(miller_set, data=rt_subset['intensity_for_norm'])
  n_refl = rt_subset.size()

  #set up binning objects
  if n_refl > 20000:
    n_refl_shells = 20
  elif n_refl > 15000:
    n_refl_shells = 15
  elif n_refl > 10000:
    n_refl_shells = 10
  else:
    logger.info(('No normalised intensity values were calculated, as an insufficient\n'
    'number of reflections were detected. All normalised intensity \n'
    'values will be set to 1 to allow use in scaling model determination. \n'))
    reflection_table['Esq'] = flex.double(reflection_table.size(), 1.0)
    return reflection_table

  d_star_sq = miller_array.d_star_sq().data()
  d_star_sq_min = flex.min(d_star_sq)
  d_star_sq_max = flex.max(d_star_sq)
  span = d_star_sq_max - d_star_sq_min

  relative_tolerance = 1e-6
  d_star_sq_max += span * relative_tolerance
  d_star_sq_min -= span * relative_tolerance
  step = (d_star_sq_max - d_star_sq_min) / n_refl_shells
  miller_array.setup_binner_d_star_sq_step(
    auto_binning=False,
    d_max=uctbx.d_star_sq_as_d(d_star_sq_max),
    d_min=uctbx.d_star_sq_as_d(d_star_sq_min),
    d_star_sq_step=step)

  normalisations = miller_array.intensity_quasi_normalisations()
  normalised_intensities = miller_array.customized_copy(
    data=(miller_array.data()/normalisations.data()))
  reflection_table['Esq'] = flex.double(reflection_table.size(), 0.0)
  reflection_table['Esq'].set_selected(
    good_refl_sel, normalised_intensities.data())
  return reflection_table

def set_wilson_outliers(reflection_table):
  """Function that takes in a reflection table with 'Esq' and 'centric_flag'
  values and sets an outlier flag depending on a cutoff for p < 1e-6."""

  centric_cutoff = 23.91
  sel1 = reflection_table['centric_flag']
  sel2 = reflection_table['Esq'] > centric_cutoff #probability <10^-6
  reflection_table.set_flags(sel1 & sel2,
    reflection_table.flags.outlier_in_scaling)

  acentric_cutoff = 13.82
  sel1 = ~reflection_table['centric_flag']
  sel2 = reflection_table['Esq'] > acentric_cutoff #probability <10^-6
  reflection_table.set_flags(sel1 & sel2,
    reflection_table.flags.outlier_in_scaling)
  msg = ('{0} reflections have been identified as outliers based on their normalised {sep}'
    'intensity values. These are reflections that have a probablity of {sep}'
    '< 10e-6 based on a Wilson distribution (E^2 > {1}, {2} for centric {sep}'
    'and acentric reflections respectively). {sep}'
    ).format(reflection_table.get_flags(
    reflection_table.flags.outlier_in_scaling).count(True),
    centric_cutoff, acentric_cutoff, sep='\n')
  logger.info(msg)
  return reflection_table



def combine_intensities(reflection_tables, experiment, Imids=None):
  """Test various combinations of prf/sum intensities to determine optimal.
  No outlier rejection is performed as it is expected that this function will
  be called after a round of outlier rejection."""

  logger.info("Performing profile/summation intensity optimisation")
  # first analyse which reflection tables have both sum and prf intensities.
  intensities, indices_to_use, indices_to_skip = _interpret_multi_tables(reflection_tables)
  if len(indices_to_skip) == len(reflection_tables):
    logger.info('No reflection tables found with both prf and sum values,'
      'no intensity combination can be performed')
    for table in reflection_tables:
      table = _set_intensities_as_sum(table)
    return reflection_tables, None

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
  rows = []
  results = {}

  if len(reflection_tables) == 1:
    #select good reflections to use
    reflections = _filter_reflections_for_combining(reflection_tables[0])
    # Calculate prescaling corrections only once
    prescaling_corrections = [calculate_prescaling_correction(reflections)]
    for Imid in Imid_list:
      Int, Var = _get_Is_from_Imidval(reflections, Imid)
      miller_set = miller.set(crystal_symmetry=experiment.crystal.get_crystal_symmetry(),
        indices=reflections['miller_index'], anomalous_flag=False)
      i_obs = miller.array(miller_set,
        data=Int * prescaling_corrections[0] /reflections['inverse_scale_factor'])
      i_obs.set_observation_type_xray_intensity()
      i_obs.set_sigmas((Var**0.5) * prescaling_corrections[0] /reflections['inverse_scale_factor'])
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

  else:
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

  st = simple_table(rows, header)
  logger.info(st.format())

  max_key = min(results, key=results.get)
  if max_key == 0:
    logger.info('Profile intensities determined to be best for scaling. \n')
  elif max_key == 1:
    logger.info('Summation intensities determined to be best for scaling. \n')
  else:
    logger.info('Combined intensities with Imid = %s determined to be best for scaling. \n' % max_key)
  # Now we know what is the best, go through all tables, combining on the best
  # value and applying the prescaling correction
  for i in indices_to_use: #choose the right values for tables with both prf and sum
    if max_key == 0:
      conv = calculate_prescaling_correction(reflection_tables[i])
      reflection_tables[i]['intensity'] = reflection_tables[i]['intensity.prf.value'] * conv
      reflection_tables[i]['variance'] = reflection_tables[i]['intensity.prf.variance'] * conv * conv
      reflection_tables[i].set_flags(reflection_tables[i]['variance'] <= 0.0,
        reflection_tables[i].flags.excluded_for_scaling)
    elif max_key == 1:
      reflection_tables[i] = _set_intensities_as_sum(reflection_tables[i])
    else: #combine
      #do partiality selection an appliction first, then use same function.
      if 'partiality' in reflection_tables[i]:
        inverse_partiality = _determine_inverse_partiality(reflection_tables[i])
        Icomb, Vcomb = _calculate_combined_raw_intensities(
          reflection_tables[i]['intensity.prf.value'],
          reflection_tables[i]['intensity.sum.value'] * inverse_partiality,
          reflection_tables[i]['intensity.prf.variance'],
          reflection_tables[i]['intensity.sum.variance'] * (inverse_partiality**2),
          max_key)
      else:
        Icomb, Vcomb = _calculate_combined_raw_intensities(
          reflection_tables[i]['intensity.prf.value'],
          reflection_tables[i]['intensity.sum.value'],
          reflection_tables[i]['intensity.prf.variance'],
          reflection_tables[i]['intensity.sum.variance'],
          max_key)
      conv = calculate_prescaling_correction(reflection_tables[i])
      reflection_tables[i]['intensity'] = Icomb * conv
      reflection_tables[i]['variance'] = Vcomb * conv * conv
      reflection_tables[i].set_flags(reflection_tables[i]['variance'] <= 0.0,
        reflection_tables[i].flags.excluded_for_scaling)
  for i in indices_to_skip:
    reflection_tables[i] = _set_intensities_as_sum(reflection_tables[i])
  return reflection_tables, results

### Helper functions for combine_intensities

def _interpret_multi_tables(reflection_tables):
  if len(reflection_tables) == 1:
    if 'intensity.sum.value' in reflection_tables[0] and 'intensity.prf.value' in reflection_tables[0]:
      return reflection_tables[0]['intensity.sum.value'], [0], []
    else:
      return [], [], [0]
  else:
    intensities = flex.double([])
    indices_to_use = []
    indices_to_skip = []
    for i, table in enumerate(reflection_tables):
      if 'intensity.sum.value' in table and 'intensity.prf.value' in table:
        intensities.extend(table['intensity.sum.value'].as_double())
        indices_to_use.append(i)
      else:
        assert 'intensity.sum.value' in table
        indices_to_skip.append(i)
    return intensities, indices_to_use, indices_to_skip

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

def _calculate_combined_raw_intensities(Iprf, Isum, Vprf, Vsum, Imid):
  """Use partiality-corrected Isum, alongside Iprf to calculate
  combined raw intensities."""
  w = 1.0/(1.0 + (Isum/Imid)**3)
  w.set_selected(Isum <= 0, 1.0)
  Icomb = (w * Iprf) + ((1.0 - w) * Isum)
  Vcomb = (w * Vprf) + ((1.0 - w) * Vsum)
  return Icomb, Vcomb

def calculate_prescaling_correction(reflection_table):
  """Calculate the multiplicative conversion factor for intensities."""
  conversion = flex.double(reflection_table.size(), 1.0)
  if 'lp' in reflection_table:
    conversion *= reflection_table['lp']
  qe = None
  if 'qe' in reflection_table:
    qe = reflection_table['qe']
  elif 'dqe' in reflection_table:
    qe = reflection_table['dqe']
  if qe:
    inverse_qe = flex.double(reflection_table.size(), 1.0)
    nonzero_qe_sel = qe > 0.0
    good_qe = qe.select(qe > 0.0)
    inverse_qe.set_selected(nonzero_qe_sel.iselection(), 1.0/good_qe)
    conversion *= inverse_qe
  return conversion
