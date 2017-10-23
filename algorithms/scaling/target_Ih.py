import numpy as np
from dials.array_family import flex

class single_Ih_table(object):
  def __init__(self, refl_table, weighting):
    #first create a minimal reflection table object
    self.Ih_table = flex.reflection_table()
    self.Ih_table['asu_miller_index'] = refl_table['asu_miller_index']
    self.Ih_table['intensity'] = refl_table['intensity']
    #bring in weights and initial scale factors
    self.weights_for_scaling = weighting.get_weights()
    self.scale_factors = refl_table['inverse_scale_factor']
    #calculate the indexing arrays
    (self.h_index_counter_array, self.h_index_cumulative_array) = self.assign_h_index()
    self.Ih_array = None
    #calculate a first estimate of Ih
    self.calc_Ih()

  def update_scale_factors(self, scalefactors):
    self.scale_factors = scalefactors

  def get_Ih_values(self):
    return self.Ih_table['Ih_values']

  def calc_Ih(self):
    '''calculate the current best estimate for I for each reflection group'''
    intensities = self.Ih_table['intensity']
    scale_factors = self.scale_factors
    scaleweights = self.weights_for_scaling
    gsq = (((scale_factors)**2) * scaleweights)
    sumgsq = flex.double(np.add.reduceat(gsq, self.h_index_cumulative_array[:-1]))
    gI = ((scale_factors * intensities) * scaleweights)
    sumgI = flex.double(np.add.reduceat(gI, self.h_index_cumulative_array[:-1]))
    sumweights = flex.double(np.add.reduceat(scaleweights, self.h_index_cumulative_array[:-1]))
    self.Ih_array = flex.double([val/ sumgsq[i] if sumweights[i] > 0.0
                                 else 0.0 for i, val in enumerate(sumgI)])
    self.Ih_table['Ih_values'] = flex.double(
      np.repeat(self.Ih_array, self.h_index_counter_array))

  def assign_h_index(self):
    '''assign an index to the sorted reflection table that
       labels each group of unique miller indices'''
    s = len(self.Ih_table)
    self.Ih_table['h_index'] = flex.int([0] * s)
    h_index_counter_array = []
    h_index = 0
    h_index_counter = 1
    for i in range(1, s):
      if (self.Ih_table['asu_miller_index'][i] ==
          self.Ih_table['asu_miller_index'][i-1]):
        self.Ih_table['h_index'][i] = h_index
        h_index_counter += 1
      else:
        h_index += 1
        self.Ih_table['h_index'][i] = h_index
        h_index_counter_array.append(h_index_counter)
        h_index_counter = 1
    h_index_counter_array.append(h_index_counter)
    '''calculate the cumulative sum after each h_index group'''
    hsum = 0
    h_index_cumulative_array = [0]
    for n in h_index_counter_array:
      hsum += n
      h_index_cumulative_array.append(hsum)
    return h_index_counter_array, h_index_cumulative_array


'''class target_Ih(object):
  def __init__(self, data_man):
    self.Ih_table = flex.reflection_table()
    self.Ih_table['asu_miller_index'] = data_man.sorted_reflections['asu_miller_index']
    self.Ih_table['intensity'] = data_man.sorted_reflections['intensity']
    self.Ih_table['variance'] = data_man.sorted_reflections['variance']

  def '''
  