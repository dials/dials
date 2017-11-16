import numpy as np
from dials.array_family import flex
from cctbx import miller, crystal
import time as time

class base_Ih_table(object):
  def __init__(self, refl_table, weights):
    #check necessary columns exists in input reflection table
    for column in ['asu_miller_index', 'intensity', 'inverse_scale_factor']:
      if not column in refl_table.keys():
        assert 0, """Attempting to create an Ih_table object from a reflection
        table with no %s column""" % column
    #first create a minimal reflection table object
    self.Ih_table = flex.reflection_table()
    self.Ih_table['asu_miller_index'] = refl_table['asu_miller_index']
    self.Ih_table['intensity'] = refl_table['intensity']
    self.Ih_table['Ih_values'] = flex.double([0.0]*len(refl_table))
    self.Ih_table['weights'] = weights
    self.Ih_table['inverse_scale_factor'] = refl_table['inverse_scale_factor']
    #calculate the indexing arrays
    (self.h_index_counter_array, self.h_index_cumulative_array) = self.assign_h_index()
    self.n_h = self.calc_nh()
    self.Ih_array = None #This may not be necessary in future but keep for now.
  
  #note: no calc_Ih method here, this must be filled in by subclasses, this is
  #necessary to allow scaling against a target Ih external to the Ih_table.

  def update_scale_factors(self, scalefactors):
    if len(scalefactors) != len(self.Ih_table['inverse_scale_factor']):
      assert 0, """attempting to set a new set of scale factors of different
      length than previous assignment: was %s, attempting %s""" % (
        len(self.Ih_table['inverse_scale_factor']), len(scalefactors))
    self.Ih_table['inverse_scale_factor'] = scalefactors

  def update_weights(self, weights):
    if len(weights) != len(self.Ih_table['weights']):
      assert 0, """attempting to set a new set of weights of different
      length than previous assignment: was %s, attempting %s""" % (
        len(self.Ih_table['weights']), len(weights))
    self.Ih_table['weights'] = weights

  def update_aimless_error_model(self, error_params):
    sigmaprime = error_params[0] * (((1.0/self.Ih_table['weights']) 
                                    + ((error_params[1] * self.Ih_table['intensity'])**2))**0.5)
    self.Ih_table['weights'] = 1.0/(sigmaprime**2)

  def set_Ih_values(self, Ih_values):
    if len(Ih_values) != len(self.Ih_table['Ih_values']):
      assert 0, """attempting to set a new set of Ih_values of different
      length than previous assignment: was %s, attempting %s""" % (
        len(self.Ih_table['Ih_values']), len(Ih_values))
    self.Ih_table['Ih_values'] = Ih_values

  def get_Ih_values(self):
    return self.Ih_table['Ih_values']

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

  def calc_nh(self):
    '''returns a vector of len(reflections) with the number of members of
    each h group'''
    n_h = flex.double([])
    for i in self.h_index_counter_array:
      n_h.extend(flex.double([i]*i))
    return n_h

class single_Ih_table(base_Ih_table):
  '''subclass of base_Ih_table to fill in the calc_Ih method. This is the default
  data structure used for scaling a single sweep.'''
  def __init__(self, refl_table, weighting):
    base_Ih_table.__init__(self, refl_table, weighting)
    self.calc_Ih() #calculate a first estimate of Ih

  def calc_Ih(self):
    '''calculate the current best estimate for I for each reflection group'''
    intensities = self.Ih_table['intensity']
    scale_factors = self.Ih_table['inverse_scale_factor']
    scaleweights = self.Ih_table['weights']
    gsq = (((scale_factors)**2) * scaleweights)
    sumgsq = flex.double(np.add.reduceat(gsq, self.h_index_cumulative_array[:-1]))
    gI = ((scale_factors * intensities) * scaleweights)
    sumgI = flex.double(np.add.reduceat(gI, self.h_index_cumulative_array[:-1]))
    sumweights = flex.double(np.add.reduceat(scaleweights, self.h_index_cumulative_array[:-1]))
    self.Ih_array = flex.double([val/ sumgsq[i] if sumweights[i] > 0.0
                                 else 0.0 for i, val in enumerate(sumgI)])
    self.Ih_table['Ih_values'] = flex.double(
      np.repeat(self.Ih_array, self.h_index_counter_array))


class target_Ih(object):
  def __init__(self, Ih_table_1, Ih_table_2, experiments):
    self.experiments = experiments
    self.Ih_table_1 = Ih_table_1
    self.Ih_table_2 = Ih_table_2
    self.Ih_table = flex.reflection_table()
    self.determine_all_unique_indices()#fill in a unique index column
    self.assign_hjoin_index()
    #self.calc_Ih()

  def determine_all_unique_indices(self):
    u_c = self.experiments.crystal.get_unit_cell().parameters()
    s_g = self.experiments.crystal.get_space_group()
    crystal_symmetry = crystal.symmetry(unit_cell=u_c, space_group=s_g)
    set1 = list(set(self.Ih_table_1.Ih_table['asu_miller_index']))
    set2 = list(set(self.Ih_table_2.Ih_table['asu_miller_index']))
    all_miller_indices = flex.miller_index(list(set(set1+set2)))
    miller_set = miller.set(crystal_symmetry=crystal_symmetry,
                            indices=all_miller_indices)
    permuted = miller_set.sort_permutation(by_value='packed_indices')
    self.Ih_table['unique_indices'] = all_miller_indices.select(permuted)

  def assign_hjoin_index(self):
    #looks at sorted miller indices and counts instances relative to the target
    miller_index_1 = list(self.Ih_table_1.Ih_table['asu_miller_index'])
    miller_index_2 = list(self.Ih_table_2.Ih_table['asu_miller_index'])
    self.h_idx_count_1 = []
    self.h_idx_count_2 = []
    for unique_index in self.Ih_table['unique_indices']:
      self.h_idx_count_1.append(miller_index_1.count(unique_index))
      self.h_idx_count_2.append(miller_index_2.count(unique_index))
    hsum_1 = 0
    hsum_2 = 0
    self.h_index_cumulative_array_1 = [0]
    self.h_index_cumulative_array_2 = [0]
    for n in self.h_idx_count_1:
      hsum_1 += n
      self.h_index_cumulative_array_1.append(hsum_1)
    for n in self.h_idx_count_2:
      hsum_2 += n
      self.h_index_cumulative_array_2.append(hsum_2)

  def calc_Ih(self):
    self.Ih_table['Ih_values'] = flex.double([0.0]*len(self.Ih_table))
    for i, _ in enumerate(self.Ih_table):
      (s1, f1) = (self.h_index_cumulative_array_1[i], self.h_index_cumulative_array_1[i+1])
      (s2, f2) = (self.h_index_cumulative_array_2[i], self.h_index_cumulative_array_2[i+1])
      weights = self.Ih_table_1.weights_for_scaling[s1:f1]
      weights.extend(self.Ih_table_2.weights_for_scaling[s2:f2])
      intensities = self.Ih_table_1.Ih_table['intensity'][s1:f1]
      intensities.extend(self.Ih_table_2.Ih_table['intensity'][s2:f2])
      scales = self.Ih_table_1.scale_factors[s1:f1]
      scales.extend(self.Ih_table_2.scale_factors[s2:f2])
      self.Ih_table['Ih_values'][i] = (flex.sum(scales*weights*intensities)
                                       / flex.sum(weights*(scales**2)))

  def return_Ih_values(self):
    Ih1_values = flex.double(np.repeat(self.Ih_table['Ih_values'], self.h_idx_count_1))
    Ih2_values = flex.double(np.repeat(self.Ih_table['Ih_values'], self.h_idx_count_2))
    return Ih1_values, Ih2_values
