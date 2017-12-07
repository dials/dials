'''
Simple functions for calculating R_meas, R_int from a Data_Manager_object
'''
from dials.array_family import flex
import numpy as np

def R_pim_meas(data_man):
  '''Calculate R_meas from a Data_Manager_object'''
  Ihl = data_man.Ih_table.Ih_table['intensity']
  gvalues = data_man.Ih_table.Ih_table['inverse_scale_factor']

  ones = flex.double([1.0] * len(Ihl))
  nh = ones * data_man.Ih_table.h_index_mat

  I_average = (((Ihl/gvalues) * data_man.Ih_table.h_index_mat)/nh)
  I_average_expanded = flex.double(np.repeat(I_average, 
    data_man.Ih_table.h_index_counter_array))

  diff = abs((Ihl/gvalues) - I_average_expanded)
  reduced_diff = diff * data_man.Ih_table.h_index_mat
 
  selection = (nh != 1.0)
  sel_reduced_diff = reduced_diff.select(selection)
  sel_nh = nh.select(selection)

  Rpim_upper = flex.sum(((1.0/(sel_nh - 1.0))**0.5) * sel_reduced_diff)
  Rmeas_upper = flex.sum(((sel_nh/(sel_nh - 1.0))**0.5) * sel_reduced_diff)
  sumIh = I_average_expanded * data_man.Ih_table.h_index_mat
  sumIh = sumIh.select(selection)
  Rpim_lower = flex.sum(sumIh)
  Rpim = Rpim_upper/Rpim_lower
  Rmeas = Rmeas_upper/Rpim_lower
  return Rpim, Rmeas