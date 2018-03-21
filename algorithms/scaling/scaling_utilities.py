from __future__ import print_function
import logging
from time import time
import copy
import cPickle as pickle
from dials.array_family import flex
from math import pi
from cctbx import miller, crystal
from Ih_table import SingleIhTable
from dxtbx.model.experiment_list import ExperimentListDumper
from dials_scratch_scaling_ext import create_sph_harm_table, calc_theta_phi,\
  rotate_vectors_about_axis

logger = logging.getLogger('dials')

def save_experiments(experiments, filename):
  ''' Save the experiments json. '''
  st = time()
  logger.info('Saving the experiments to %s' % filename)
  dump = ExperimentListDumper(experiments)
  with open(filename, "w") as outfile:
    outfile.write(dump.as_json(split=True))
  logger.info('Time taken: %g' % (time() - st))

def save_reflections(scaler, filename):
  '''Save the scaled reflections.'''
  #First remove any additional columns added during scaling that are not needed.
  scaler.clean_reflection_table() 
  st = time()
  logger.info('Saving the scaled reflections to %s' % filename)
  scaler.reflection_table.as_pickle(filename)
  logger.info('Time taken: %g' % (time() - st))

def parse_multiple_datasets(reflections):
  'method to parse multiple datasets from single reflection tables, selecting on id'
  single_reflection_tables = []
  for refl_table in reflections:
    dataset_ids = set(refl_table['id']).difference(set([-1]))
    n_datasets = len(dataset_ids)
    if n_datasets > 1:
      logger.info(('\nDetected existence of a multi-dataset scaled reflection table {sep}'
        'containing {0} datasets. {sep}').format(n_datasets, sep='\n'))
      for dataset_id in dataset_ids:
        single_refl_table = refl_table.select(refl_table['id'] == dataset_id)
        single_reflection_tables.append(single_refl_table)
    else:
      single_reflection_tables.append(refl_table)
  return single_reflection_tables

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
  (ux, uy, uz) = exp_rot_axis[0][0], exp_rot_axis[0][1], exp_rot_axis[0][2]
  cross_prod_uz = flex.vec3_double([(uy, -1.0*ux, 0.0)])
  from math import acos
  angle_between_u_z = -1.0*acos(uz/((ux**2 + uy**2 + uz**2)**0.5))
  phi = flex.double([angle_between_u_z]*len(vectors))
  new_vectors = rotate_vectors_about_axis(cross_prod_uz, vectors, phi)
  return flex.vec3_double(new_vectors)

def sph_harm_table(reflection_table, experiments, lmax):
  """Calculate the spherical harmonic table for a spherical
    harmonic absorption correction."""
  from scitbx import sparse # Leave - needed to allow return of sph_harm_table
  reflection_table = calc_crystal_frame_vectors(reflection_table, experiments)
  theta_phi = calc_theta_phi(reflection_table['s0c'])
  theta_phi_2 = calc_theta_phi(reflection_table['s1c'])
  sph_h_t = create_sph_harm_table(theta_phi, theta_phi_2, lmax)
  return sph_h_t

def reject_outliers(reflection_table, zmax):
  '''simple, quick, outlier rejection based on normalised deviations
  (similar to aimless)'''
  Ih_table = SingleIhTable(reflection_table)
  I = Ih_table.intensities
  g = Ih_table.inverse_scale_factors
  w = Ih_table.weights
  wgIsum = ((w * g * I) * Ih_table.h_index_matrix) * Ih_table.h_expand_matrix
  wg2sum = ((w * g * g) * Ih_table.h_index_matrix) * Ih_table.h_expand_matrix
  norm_dev = (I - (g * wgIsum/wg2sum))/(((1.0/w)+((g/wg2sum)**2))**0.5)
  z_score = (norm_dev**2)**0.5
  outliers_sel = z_score > zmax
  select_isel = Ih_table._nonzero_weights.iselection()
  outliers_in_overall = select_isel.select(outliers_sel)
  outlier_mask = flex.bool([False]*len(reflection_table))
  outlier_mask.set_selected(outliers_in_overall, True)
  reflection_table.set_flags(outlier_mask,
    reflection_table.flags.outlier_in_scaling)
  return reflection_table


def calc_normE2(reflection_table, experiments):
  '''calculate normalised intensity values for centric and acentric reflections'''
  msg = ('Calculating normalised intensity values for outlier rejection. {sep}'
    'Negative intensities are set to zero for the purpose of calculating {sep}'
    'mean intensity values for resolution bins. This is to avoid spuriously {sep}'
    'high E^2 values due to a mean close to zero and should only affect {sep}'
    'the E^2 values of the highest resolution bins. {sep}'
    ).format(sep='\n')
  logger.info(msg)
  #print(len(reflection_table))
  sel = reflection_table.get_flags(reflection_table.flags.bad_for_scaling,
    all=False)
  good_refl = ~sel
  scaling_subset = reflection_table.select(good_refl)
  # Scaling subset is data that has not been flagged as bad or excluded
  u_c = experiments.crystal.get_unit_cell().parameters()
  s_g = experiments.crystal.get_space_group()
  crystal_symmetry = crystal.symmetry(unit_cell=u_c, space_group=s_g)
  miller_set = miller.set(crystal_symmetry=crystal_symmetry,
                          indices=scaling_subset['asu_miller_index'])
  scaling_subset['resolution'] = 1.0/scaling_subset['d']**2
  #handle negative reflections to minimise effect on mean I values.
  scaling_subset['intensity_for_norm'] = copy.deepcopy(scaling_subset['intensity'])
  sel = scaling_subset['intensity'] < 0.0
  scaling_subset['intensity_for_norm'].set_selected(sel, 0.0)
  miller_array = miller.array(miller_set, data=scaling_subset['intensity_for_norm'])
  #set up binning objects
  scaling_subset['centric_flag'] = miller_array.centric_flags().data()
  n_centrics = scaling_subset['centric_flag'].count(True)
  n_acentrics = scaling_subset['centric_flag'].count(False)

  if n_acentrics > 20000 or n_centrics > 20000:
    n_refl_shells = 20
  elif n_acentrics > 15000 or n_centrics > 15000:
    n_refl_shells = 15
  elif n_acentrics < 10000:
    reflection_table['Esq'] = flex.double([1.0]*len(reflection_table))
    del reflection_table['intensity_for_norm']
    msg = ('No normalised intensity values were calculated, as an {sep}'
    'insufficient number of reflections were detected. {sep}'
    ).format(sep='\n')
    logger.info(msg)
    return reflection_table
  else:
    n_refl_shells = 10

  #calculate normalised intensities: first calculate bin averages
  step = ((max(scaling_subset['resolution']) - min(scaling_subset['resolution'])
           + 1e-8) / n_refl_shells)
  if n_centrics:
    centrics_array = miller_array.select_centric()
    centric_binner = centrics_array.setup_binner_d_star_sq_step(d_star_sq_step=step)
    mean_centric_values = centrics_array.mean(use_binning=centric_binner)
    mean_centric_values = mean_centric_values.data[1:-1]
    centric_bin_limits = centric_binner.limits()
  #  print(list(centric_bin_limits))

  acentrics_array = miller_array.select_acentric()
  acentric_binner = acentrics_array.setup_binner_d_star_sq_step(d_star_sq_step=step)
  mean_acentric_values = acentrics_array.mean(use_binning=acentric_binner)
  mean_acentric_values = mean_acentric_values.data[1:-1]
  acentric_bin_limits = acentric_binner.limits()
  #print(list(acentric_bin_limits))
  #now calculate normalised intensity values for full reflection table
  miller_set = miller.set(crystal_symmetry=crystal_symmetry,
                          indices=reflection_table['asu_miller_index'])

  reflection_table['Esq'] = flex.double([0.0]*len(reflection_table))
  miller_array = miller.array(miller_set)
  reflection_table['centric_flag'] = miller_array.centric_flags().data()
  d0_sel = reflection_table['d'] == 0.0
  reflection_table['d'].set_selected(d0_sel, 1.0)  #set for now, then set back to zero later
  reflection_table['resolution'] = 1.0/reflection_table['d']**2

  if n_centrics:
    for i in range(0, len(centric_bin_limits)-1):
      sel1 = reflection_table['centric_flag'] == True
      sel2 = reflection_table['resolution'] > centric_bin_limits[i]
      sel3 = reflection_table['resolution'] <= centric_bin_limits[i+1]
      sel = sel1 & sel2 & sel3
      intensities = reflection_table['intensity'].select(sel)
      reflection_table['Esq'].set_selected(sel, intensities/ mean_centric_values[i])
  for i in range(0, len(acentric_bin_limits)-1):
    sel1 = reflection_table['centric_flag'] == False
    sel2 = reflection_table['resolution'] > acentric_bin_limits[i]
    sel3 = reflection_table['resolution'] <= acentric_bin_limits[i+1]
    sel = sel1 & sel2 & sel3
    intensities = reflection_table['intensity'].select(sel)
    reflection_table['Esq'].set_selected(sel, intensities/ mean_acentric_values[i])
  del reflection_table['intensity_for_norm']
  reflection_table['d'].set_selected(d0_sel, 0.0)
  reflection_table['Esq'].set_selected(d0_sel, 0.0)
  msg = ('The number of centric & acentric reflections is {0} & {1}, {sep}'
    '{2} resolution bins were used for the E^2 calculation. {sep}'
    ).format(n_centrics, n_acentrics, n_refl_shells, sep='\n')
  logger.info(msg)
  return reflection_table

def calculate_wilson_outliers(reflection_table):
  """function that takes in a reflection table and experiments object and
  looks at the wilson distribution of intensities in reflection shells to
  look for the presence of outliers with high intensities. Returns a bool
  flex array indicating any outliers."""

  centric_cutoff = 23.91
  sel1 = reflection_table['centric_flag'] == True
  sel2 = reflection_table['Esq'] > centric_cutoff #probability <10^-6
  reflection_table.set_flags(sel1 & sel2, reflection_table.flags.outlier_in_scaling)

  acentric_cutoff = 13.82
  sel1 = reflection_table['centric_flag'] == False
  sel2 = reflection_table['Esq'] > acentric_cutoff #probability <10^-6
  reflection_table.set_flags(sel1 & sel2, reflection_table.flags.outlier_in_scaling)
  msg = ('{0} reflections have been identified as outliers based on their normalised {sep}'
    'intensity values. These are reflections that have a probablity of {sep}'
    '< 10e-6 based on a Wilson distribution (E^2 > {1}, {2} for centric {sep}'
    'and acentric reflections respectively). {sep}'
    ).format(reflection_table.get_flags(
    reflection_table.flags.outlier_in_scaling).count(True),
    centric_cutoff,acentric_cutoff, sep='\n')
  logger.info(msg)
  return reflection_table

'''def R_pim_meas(scaler):
  #Calculate R_pim, R_meas from a scaler
  Ihl = scaler.Ih_table.intensities
  gvalues = scaler.Ih_table.inverse_scale_factors

  ones = flex.double([1.0] * len(Ihl))
  nh = ones * scaler.Ih_table.h_index_matrix

  I_average = (((Ihl/gvalues) * scaler.Ih_table.h_index_matrix)/nh)
  I_average_expanded = flex.double(np.repeat(I_average,
    scaler.Ih_table.h_index_counter_array))

  diff = abs((Ihl/gvalues) - I_average_expanded)
  reduced_diff = diff * scaler.Ih_table.h_index_matrix

  selection = (nh != 1.0)
  sel_reduced_diff = reduced_diff.select(selection)
  sel_nh = nh.select(selection)

  Rpim_upper = flex.sum(((1.0/(sel_nh - 1.0))**0.5) * sel_reduced_diff)
  Rmeas_upper = flex.sum(((sel_nh/(sel_nh - 1.0))**0.5) * sel_reduced_diff)
  sumIh = I_average_expanded * scaler.Ih_table.h_index_matrix
  sumIh = sumIh.select(selection)
  Rpim_lower = flex.sum(sumIh)
  Rpim = Rpim_upper/Rpim_lower
  Rmeas = Rmeas_upper/Rpim_lower
  return Rpim, Rmeas'''
