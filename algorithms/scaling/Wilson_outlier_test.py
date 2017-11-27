import copy
from dials.array_family import flex
from cctbx import miller, crystal
import numpy as np
from target_function import *
from basis_functions import *
import scale_factor as SF
from reflection_weighting import *
from data_quality_assessment import R_meas, R_pim
import matplotlib.pyplot as plt

def calc_normE2(reflection_table, experiments):
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
  reflection_table['centric_flag'] = miller_array.centric_flags().data()
  #set up binning objects
  n_centrics = reflection_table['centric_flag'].count(True)
  n_acentrics = reflection_table['centric_flag'].count(False)
  print """Number of centric/acentric reflections is
           %s, %s""" % (n_centrics, n_acentrics)
  if n_acentrics > 20000 or n_centrics > 20000:
    n_refl_shells = 20
  elif n_acentrics > 15000 or n_centrics > 15000:
    n_refl_shells = 15
  else:
    n_refl_shells = 10
  #calculate normalised intensities
  #first bin by resolution and determine the average values of each bin

  step = ((max(reflection_table['resolution']) - min(reflection_table['resolution'])
           + 1e-8) / n_refl_shells)
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

  reflection_table['Esq'] = flex.double([0.0]*len(reflection_table))
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
  return reflection_table

'''def calculate_normalised_intensities(reflection_table):
  #want each reflection 'shell' to contain at least 500 reflections.'''
  


def calculate_wilson_outliers(reflection_table):
  '''function that takes in a reflection table and experiments object and
  looks at the wilson distribution of intensities in reflection shells to
  look for the presence of outliers with high intensities. Returns a bool
  flex array indicating any outliers.'''
  reflection_table['wilson_outlier_flag'] = flex.bool([False] * len(reflection_table))

  sel1 = reflection_table['centric_flag'] == True
  sel2 = reflection_table['Esq'] > 23.91 #probability <10^-6
  reflection_table['wilson_outlier_flag'].set_selected(sel1 & sel2, True)

  sel1 = reflection_table['centric_flag'] == False
  sel2 = reflection_table['Esq'] > 13.82 #probability <10^-6
  reflection_table['wilson_outlier_flag'].set_selected(sel1 & sel2, True)

  print "found %s outliers from analysis of Wilson statistics" % (
    reflection_table['wilson_outlier_flag'].count(True))
  return reflection_table['wilson_outlier_flag']

def select_variables_in_range(variable_array, lower_limit, upper_limit):
  '''return boolean selection of a given variable range'''
  sel = flex.bool()
  for variable in variable_array:
    if lower_limit < variable and variable <= upper_limit:
      sel.append(True)
    else:
      sel.append(False)
  return sel