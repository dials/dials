'''
Simple functions for calculating R_meas, R_int from a Data_Manager_object
'''

def R_meas(Data_Manager_object):
  '''Calculate R_meas from a Data_Manager_object'''
  Ihl = Data_Manager_object.sorted_reflections['intensity']
  Ih = Data_Manager_object.sorted_reflections['Ih_values']

  gvalues = Data_Manager_object.sorted_reflections['inverse_scale_factor']

  Rmeas_upper = 0.0
  Rmeas_lower = 0.0
  Ih_average = []
  for h in range(len(Data_Manager_object.h_index_counter_array)):
    a1 = 0.0
    lsum = Data_Manager_object.h_index_counter_array[h]
    if lsum > 1:
      for i in range(lsum):
        indexer = i + Data_Manager_object.h_index_cumulative_array[h]
        a1 += (Ihl[indexer]/ (gvalues[indexer]))
      average = a1/lsum
      for i in range(lsum):
        Ih_average.append(average)
    else:
      Ih_average.append(0.0)
  for h in range(len(Data_Manager_object.h_index_counter_array)):
    a1 = 0.0
    lsum = Data_Manager_object.h_index_counter_array[h]
    if lsum > 1:
      for i in range(lsum):
        indexer = i + Data_Manager_object.h_index_cumulative_array[h]
        a1 += abs((Ihl[indexer] / (gvalues[indexer])) - Ih_average[indexer])
        Rmeas_lower += (Ih_average[indexer])
      Rmeas_upper += (((float(lsum) / (float(lsum) - 1.0))**0.5) * a1)
  Rmeas = Rmeas_upper / Rmeas_lower
  return Rmeas

def R_pim(Data_Manager_object):
  '''Calculate R_pim from a Data_Manager_object'''
  Ihl = Data_Manager_object.sorted_reflections['intensity']
  Ih = Data_Manager_object.sorted_reflections['Ih_values']
  gvalues = Data_Manager_object.sorted_reflections['inverse_scale_factor']
  Rpim_upper = 0.0
  Rpim_lower = 0.0
  Ih_average=[]
  #calculate the average Ih for each group of reflections
  for h in range(len(Data_Manager_object.h_index_counter_array)):
    a1 = 0.0
    lsum = Data_Manager_object.h_index_counter_array[h]
    if lsum > 1:
      for i in range(lsum):
        indexer = i + Data_Manager_object.h_index_cumulative_array[h]
        a1 += (Ihl[indexer]/gvalues[indexer])
      average = a1/lsum
      for i in range(lsum):
        Ih_average.append(average)
    else:
      Ih_average.append(0.0)
  for h in range(len(Data_Manager_object.h_index_counter_array)):
    a1 = 0.0
    lsum = Data_Manager_object.h_index_counter_array[h]
    if lsum > 1:
      for i in range(lsum):
        indexer = i + Data_Manager_object.h_index_cumulative_array[h]
        a1 += abs((Ihl[indexer] / (gvalues[indexer])) - Ih_average[indexer])
        Rpim_lower += (Ih_average[indexer])
      Rpim_upper += (((1.0 / (float(lsum) - 1.0))**0.5) * a1)
  Rpim = Rpim_upper / Rpim_lower
  return Rpim
