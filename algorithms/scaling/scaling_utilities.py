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

def apply_prescaling_correction(reflection_table, conv):
  reflection_table['intensity'] *= conv
  reflection_table['variance'] *= conv * conv
  return reflection_table

def combine_intensities(reflection_tables, experiment, Imids=None):
  """Test various combinations of prf/sum intensities to determine optimal.
  No outlier rejection is performed as it is expected that this function will
  be called after a round of outlier rejection."""

  # first analyse which reflection tables have both sum and prf intensities.
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
  if len(indices_to_skip) == len(reflection_tables):
    logger.info('No reflection tables found with both prf and sum values,'
      'no intensity combination can be performed')
    for i, table in enumerate(reflection_tables):
      reflection_tables[i]['intensity'] = reflection_tables[i]['intensity.sum.value']
    return reflection_tables, None

  if Imids:
    Imid_list = Imids
  else:
    avg = flex.mean(intensities)
    sorted_intensities = flex.sorted(intensities, reverse=True)
    Imid = sorted_intensities[0]/10.0
    Imid_list = [0, 1, avg, Imid]
    while (Imid > avg):
      Imid /= 10.0
      Imid_list.append(Imid)

  # Calculate prescaling corrections only once
  prescaling_corrections = [None]*len(reflection_tables)
  for i in indices_to_use:
    reflections = reflection_tables[i]
    reflections = reflections.select(~reflections.get_flags(
      reflections.flags.bad_for_scaling, all=False))
    reflections = reflections.select(reflections['intensity.prf.variance'] > 0)
    reflections = reflections.select(reflections['intensity.sum.variance'] > 0)
    if 'partiality' in reflections:
      reflections = reflections.select(reflections['partiality'] > 0)
    prescaling_corrections[i] = calculate_prescaling_correction(reflections)

  header = ['Combination', 'CC1/2', 'Rmeas']
  rows = []
  results = {}

  for Is in Imid_list:
    combined_intensities = flex.double([])
    combined_variances = flex.double([])
    combined_scales = flex.double([])
    combined_indices = flex.miller_index([])
    #calculate combined intensities
    for i in indices_to_use:
      reflections = reflection_tables[i]
      #do filtering
      reflections = reflections.select(~reflections.get_flags(
        reflections.flags.bad_for_scaling, all=False))
      reflections = reflections.select(reflections['intensity.prf.variance'] > 0)
      reflections = reflections.select(reflections['intensity.sum.variance'] > 0)
      if 'partiality' in reflections:
        reflections = reflections.select(reflections['partiality'] > 0)
      # do calculation and add to data
      if Is == 0: #special value to trigger prf
        reflections['intensity'] = reflections['intensity.prf.value']
        reflections['variance'] = reflections['intensity.prf.variance']
      elif Is == 1: #special value to trigger sum
        if 'partiality' in reflections:
          reflections['intensity'] = reflections['intensity.sum.value']/reflections['partiality']
          reflections['variance'] = reflections['intensity.sum.variance']/(reflections['partiality']**2)
        else:
          reflections['intensity'] = reflections['intensity.sum.value']
          reflections['variance'] = reflections['intensity.sum.variance']
      else:
        reflections = calculate_combined_raw_intensities(reflections, Is)
      reflections = apply_prescaling_correction(reflections, prescaling_corrections[i])
      combined_intensities.extend(reflections['intensity'])
      combined_variances.extend(reflections['variance'])
      combined_scales.extend(reflections['inverse_scale_factor'])
      combined_indices.extend(reflections['miller_index'])
    # now can calculate combined statistics
    miller_set = miller.set(crystal_symmetry=experiment.crystal.get_crystal_symmetry(),
      indices=combined_indices, anomalous_flag=False)
    i_obs = miller.array(miller_set, data=combined_intensities/combined_scales)
    i_obs.set_observation_type_xray_intensity()
    i_obs.set_sigmas((combined_variances**0.5)/combined_scales)
    n_bins = min(20, int(combined_intensities.size()/100)+1)
    try:
      res = iotbx.merging_statistics.dataset_statistics(i_obs=i_obs, n_bins=n_bins,
        anomalous=False, sigma_filtering=None, use_internal_variance=False,
        eliminate_sys_absent=False)
    except RuntimeError:
      raise DialsMergingStatisticsError("Unable to merge for intensity combination")
    # record the results
    results[Is] = res.overall.r_meas
    if Is == 0:
      res_str = 'prf only'
    elif Is == 1:
      res_str = 'sum only'
    else:
      res_str = 'Imid = '+str(round(Is, 2))
    rows.append([res_str, str(round(res.overall.cc_one_half, 5)),
      str(round(res.overall.r_meas, 5))])

  st = simple_table(rows, header)
  logger.info(st.format())

  max_key = min(results, key=results.get)
  if max_key == 0:
    logger.info('prf intensities determined to be best for scaling. \n')
  elif max_key == 1:
    logger.info('sum intensities determined to be best for scaling. \n')
  else:
    logger.info('Combined intensities with Imid = %s determined to be best for scaling. \n' % max_key)
  # Now we know what is the best, go through all tables, combining on the best
  # value and applying the prescaling correction
  for i in indices_to_use: #choose the right values for tables with both prf and sum
    if max_key == 0:
      reflection_tables[i]['intensity'] = reflection_tables[i]['intensity.prf.value']
      reflection_tables[i]['variance'] = reflection_tables[i]['intensity.prf.variance']
    elif max_key == 1:
      inverse_partiality = flex.double(reflection_tables[i].size(), 1.0)
      if 'partiality' in reflection_tables[i]:
        nonzero_partiality_sel = reflection_tables[i]['partiality'] > 0.0
        good_refl = reflection_tables[i].select(nonzero_partiality_sel)
        inverse_partiality.set_selected(nonzero_partiality_sel.iselection(),
          1.0/good_refl['partiality'])
      reflection_tables[i]['intensity'] = \
        reflection_tables[i]['intensity.sum.value'] * inverse_partiality
      reflection_tables[i]['variance'] = \
        reflection_tables[i]['intensity.sum.variance'] * (inverse_partiality**2)
      reflection_tables[i].set_flags(reflection_tables[i]['variance'] <= 0.0,
        reflection_tables[i].flags.excluded_for_scaling)
    else:
      reflection_tables[i] = calculate_combined_raw_intensities(
        reflection_tables[i], max_key)
      reflection_tables[i].set_flags(reflection_tables[i]['variance'] <= 0.0,
        reflection_tables[i].flags.excluded_for_scaling)
    conv = calculate_prescaling_correction(reflection_tables[i])
    reflection_tables[i] = apply_prescaling_correction(reflection_tables[i], conv)
  for i in indices_to_skip:
    inverse_partiality = flex.double(reflection_tables[i].size(), 1.0)
    if 'partiality' in reflection_tables[i]:
      nonzero_partiality_sel = reflection_tables[i]['partiality'] > 0.0
      good_refl = reflection_tables[i].select(nonzero_partiality_sel)
      inverse_partiality.set_selected(nonzero_partiality_sel.iselection(),
        1.0/good_refl['partiality'])
    reflection_tables[i]['intensity'] = \
      reflection_tables[i]['intensity.sum.value'] * inverse_partiality
    reflection_tables[i]['variance'] = \
      reflection_tables[i]['intensity.sum.variance'] * (inverse_partiality**2)
    reflection_tables[i].set_flags(reflection_tables[i]['variance'] <= 0.0,
      reflection_tables[i].flags.excluded_for_scaling)
    conv = calculate_prescaling_correction(reflection_tables[i])
    reflection_tables[i] = apply_prescaling_correction(reflection_tables[i], conv)
  return reflection_tables, results

def calculate_combined_raw_intensities(reflection_table, Imid):
  #if Isum >> Imid, W > 0, intensity is sum intenisty
  w = 1.0/(1.0 + (reflection_table['intensity.sum.value']/Imid)**3)
  w.set_selected(reflection_table['intensity.sum.value'] <= 0, 1.0)
  inverse_partiality = flex.double(reflection_table.size(), 1.0)
  if 'partiality' in reflection_table:
    nonzero_partiality_sel = reflection_table['partiality'] > 0.0
    good_refl = reflection_table.select(nonzero_partiality_sel)
    inverse_partiality.set_selected(nonzero_partiality_sel.iselection(),
      1.0/good_refl['partiality'])
  reflection_table['intensity'] = (w * reflection_table['intensity.prf.value']) + \
    ((1.0 - w) * reflection_table['intensity.sum.value'] * inverse_partiality)
  reflection_table['variance'] = (w * reflection_table['intensity.prf.variance']) + \
    ((1.0 - w) * reflection_table['intensity.sum.variance']  * inverse_partiality)
  return reflection_table

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
