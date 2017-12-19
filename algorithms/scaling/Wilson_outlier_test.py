from __future__ import print_function
import copy
from dials.array_family import flex
from cctbx import miller, crystal
from reflection_weighting import Weighting

import logging
logger = logging.getLogger('dials.scale')

def calc_normE2(reflection_table, experiments):
  '''calculate normalised intensity values for centric and acentric reflections'''
  msg = ('Calculating normalised intensity values. {sep}'
    'Negative intensities are set to zero for the purpose of calculating {sep}'
    'mean intensity values for resolution bins. This is to avoid spuriously {sep}'
    'high E^2 values due to a mean close to zero and should only affect {sep}'
    'the E^2 values of the highest resolution bins. {sep}'
    ).format(sep='\n')
  logger.info(msg)
  u_c = experiments.crystal.get_unit_cell().parameters()
  s_g = experiments.crystal.get_space_group()
  crystal_symmetry = crystal.symmetry(unit_cell=u_c, space_group=s_g)
  miller_set = miller.set(crystal_symmetry=crystal_symmetry,
                          indices=reflection_table['asu_miller_index'])
  reflection_table = reflection_table.select(reflection_table['d'] > 0.0)
  reflection_table['resolution'] = 1.0/reflection_table['d']**2
  #handle negative reflections to minimise effect on mean I values.
  reflection_table['intensity_for_norm'] = copy.deepcopy(reflection_table['intensity'])
  sel = reflection_table['intensity'] < 0.0
  reflection_table['intensity_for_norm'].set_selected(sel, 0.0)
  miller_array = miller.array(miller_set, data=reflection_table['intensity_for_norm'])
  #set up binning objects
  reflection_table['centric_flag'] = miller_array.centric_flags().data()
  n_centrics = reflection_table['centric_flag'].count(True)
  n_acentrics = reflection_table['centric_flag'].count(False)
  
  if n_acentrics > 20000 or n_centrics > 20000:
    n_refl_shells = 20
  elif n_acentrics > 15000 or n_centrics > 15000:
    n_refl_shells = 15
  elif n_acentrics < 100:
    reflection_table['Esq'] = flex.double([1.0]*len(reflection_table))
    del reflection_table['intensity_for_norm']
    msg = ('No normalised intensity values were calculated, {sep}'
    'as an insufficient number of reflections were detected. {sep}'
    ).format(sep='\n')
    logger.info(msg)
    return reflection_table
  else:
    n_refl_shells = 10
  
  #calculate normalised intensities: first calculate bin averages
  step = ((max(reflection_table['resolution']) - min(reflection_table['resolution'])
           + 1e-8) / n_refl_shells)
  if n_centrics:
    centrics_array = miller_array.select_centric()
    centric_binner = centrics_array.setup_binner_d_star_sq_step(d_star_sq_step=step)
    mean_centric_values = centrics_array.mean(use_binning=centric_binner)
    mean_centric_values = mean_centric_values.data[1:-1]
    centric_bin_limits = centric_binner.limits()
      
  acentrics_array = miller_array.select_acentric()
  acentric_binner = acentrics_array.setup_binner_d_star_sq_step(d_star_sq_step=step)
  mean_acentric_values = acentrics_array.mean(use_binning=acentric_binner)
  mean_acentric_values = mean_acentric_values.data[1:-1]
  acentric_bin_limits = acentric_binner.limits()
  #now calculate normalised intensity values
  reflection_table['Esq'] = flex.double([0.0]*len(reflection_table))
  if n_centrics:
    for i in range(0,len(centric_bin_limits)-1):
      sel1 = reflection_table['centric_flag'] == True
      sel2 = reflection_table['resolution'] > centric_bin_limits[i]
      sel3 = reflection_table['resolution'] <= centric_bin_limits[i+1]
      sel = sel1 & sel2 & sel3
      intensities = reflection_table['intensity'].select(sel)
      reflection_table['Esq'].set_selected(sel, intensities/ mean_centric_values[i])
  for i in range(0,len(acentric_bin_limits)-1):
    sel1 = reflection_table['centric_flag'] == False
    sel2 = reflection_table['resolution'] > acentric_bin_limits[i]
    sel3 = reflection_table['resolution'] <= acentric_bin_limits[i+1]
    sel = sel1 & sel2 & sel3
    intensities = reflection_table['intensity'].select(sel)
    reflection_table['Esq'].set_selected(sel, intensities/ mean_acentric_values[i])
  del reflection_table['intensity_for_norm']
  msg = ('Calculated normalised intensity values. {sep}'
    'The number of centric & acentric reflections is {0} & {1}. {sep}'
    'Intensities were binned into {2} resolution bins. {sep}'
    "Normalised intensities were added to the reflection table as 'Esq'. {sep}"
    ).format(n_centrics, n_acentrics, n_refl_shells, sep='\n')
  logger.info(msg)
  return reflection_table

def calculate_wilson_outliers(reflection_table):
  '''function that takes in a reflection table and experiments object and
  looks at the wilson distribution of intensities in reflection shells to
  look for the presence of outliers with high intensities. Returns a bool
  flex array indicating any outliers.'''
  reflection_table['wilson_outlier_flag'] = flex.bool([False] * len(reflection_table))

  centric_cutoff = 23.91
  sel1 = reflection_table['centric_flag'] == True
  sel2 = reflection_table['Esq'] > centric_cutoff #probability <10^-6
  reflection_table['wilson_outlier_flag'].set_selected(sel1 & sel2, True)

  acentric_cutoff = 13.82
  sel1 = reflection_table['centric_flag'] == False
  sel2 = reflection_table['Esq'] > acentric_cutoff #probability <10^-6
  reflection_table['wilson_outlier_flag'].set_selected(sel1 & sel2, True)
  msg = ('{0} reflections were rejected as outliers based on their normalised {sep}'
    'intensity values. These are reflections that have a probablity of {sep}'
    '< 10e-6 based on a Wilson distribution (E^2 > {1}, {2} for centric {sep}'
    'and acentric reflections respectively). {sep}'
    ).format(reflection_table['wilson_outlier_flag'].count(True), centric_cutoff,
    acentric_cutoff, sep='\n')
  logger.info(msg)
  return reflection_table['wilson_outlier_flag']
