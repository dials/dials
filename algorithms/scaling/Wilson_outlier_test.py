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

def calculate_wilson_outliers(reflection_table, experiments):
    '''function that takes in a reflection table and experiments object and
    looks at the wilson distribution of intensities in reflection shells to
    look for the presence of outliers with high intensities. Returns a bool
    flex array indicating any outliers.'''
    #first create a miller_array object to get the centric flags.
    u_c = experiments.crystal.get_unit_cell().parameters()
    s_g = experiments.crystal.get_space_group()
    crystal_symmetry = crystal.symmetry(unit_cell=u_c, space_group=s_g)
    miller_set = miller.set(crystal_symmetry=crystal_symmetry,
                            indices=reflection_table['asu_miller_index'])
    miller_array = miller.array(miller_set, data=reflection_table['intensity'])
    centrics = miller_array.centric_flags()
    #neater way to do the next step -
    #problem is that one needs to extract just flags from centrics_list
    centrics_list = list(centrics)
    reflection_table['centric_flag'] = flex.bool([False] * len(reflection_table))
    reflection_table['wilson_outlier_flag'] = flex.bool([False] * len(reflection_table))
    for i in range(centrics.size()):
      reflection_table['centric_flag'][i] = centrics_list[i][1]

    #want each reflection 'shell' to contain at least 500 reflections.
    n_centrics = reflection_table['centric_flag'].count(True)
    n_acentrics = reflection_table['centric_flag'].count(False)
    #print """Number of centric/acentric reflections is
    #         %s, %s""" % (n_centrics, n_acentrics)
    if n_acentrics > 20000 or n_centrics > 20000:
      n_refl_shells = 20
    elif n_acentrics > 15000 or n_centrics > 15000:
      n_refl_shells = 15
    else:
      n_refl_shells = 10

    #calculate a resolution_bin_index for each reflections
    resmax = (1.0 / (min(reflection_table['d'])**2))
    resmin = (1.0 / (max(reflection_table['d'])**2))
    resolution_bins = ((flex.double(range(0, n_refl_shells + 1))
                        * ((resmax - resmin) / n_refl_shells))
                       + flex.double([resmin] * (n_refl_shells + 1)))
    d_bins = (1.0 /(resolution_bins[::-1]**0.5))
    #add a small tolerance to make sure no rounding errors cause extreme
    #data to not be selected
    d_bins[0] = d_bins[0] - 0.00001
    d_bins[-1] = d_bins[-1] + 0.00001
    resbin_index = flex.int([-1] * len(reflection_table['d']))
    for i in range(n_refl_shells):
      selection = select_variables_in_range(
        reflection_table['d'], d_bins[i], d_bins[i+1])
      resbin_index.set_selected(selection, i)

    #now take data in res bins - first do acentric reflections
    data_binned_into_resolution = []
    centrics_average_list = []
    acentrics_average_list = []
    for i in range(n_refl_shells):
      sel1 = resbin_index == i
      sel2 = reflection_table['centric_flag'] == False
      data_binned_into_resolution.append(reflection_table.select(sel1 & sel2))
    for i, table in enumerate(data_binned_into_resolution):
      #only do outlier test if more than 250 reflections in each shell
      if len(table) > 250:
      #handle negative intensities by setting them to zero
        sel = table['intensity'] < 0.0
        table['intensity'].set_selected(sel, 0.0)
        average = flex.mean(table['intensity'])
      else:
        average = 0.0
      acentrics_average_list.append(average)
    data_binned_into_resolution = []
    for i in range(n_refl_shells):
      sel1 = resbin_index == i
      sel2 = reflection_table['centric_flag'] == True
      data_binned_into_resolution.append(reflection_table.select(sel1 & sel2))
    for i, table in enumerate(data_binned_into_resolution):
      #only do outlier test if more than 250 reflections in each shell
      if len(table) > 250:
        #handle negative intensities by setting them to zero
        sel = table['intensity'] < 0.0
        table['intensity'].set_selected(sel, 0.0)
        average = flex.mean(table['intensity'])
      else:
        average = 0.0
      centrics_average_list.append(average)
    counter = 0
    for i, bin_index in enumerate(resbin_index):
      if reflection_table['centric_flag'][i] == False:
        if acentrics_average_list[bin_index] > 0.0:
          if (reflection_table['intensity'][i]
              / acentrics_average_list[bin_index]) > 13.82: #probability <10^-6
            reflection_table['wilson_outlier_flag'][i] = True
            counter += 1
      else:
        if centrics_average_list[bin_index] > 0.0:
          if (reflection_table['intensity'][i]
              / centrics_average_list[bin_index]) > 23.91: #probability <10^-6
            reflection_table['wilson_outlier_flag'][i] = True
            counter += 1
    print "found %s outliers from analysis of Wilson statistics" % (counter)
    return reflection_table['wilson_outlier_flag']

def select_variables_in_range(variable_array, lower_limit, upper_limit):
  '''return boolean selection of a given variable range'''
  sel = flex.bool()
  for variable in variable_array:
    if lower_limit < variable <= upper_limit:
      sel.append(True)
    else:
      sel.append(False)
  return sel