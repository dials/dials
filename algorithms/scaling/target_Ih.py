import numpy as np
from dials.array_family import flex
from cctbx import miller, crystal
import time as time
from scitbx import sparse

class base_Ih_table(object):
  def __init__(self, refl_table, weights):
    self.Ih_table = self.create_Ih_table(refl_table, weights)
    (self.h_index_counter_array, self.h_index_cumulative_array) = self.assign_h_index()
    self.h_index_mat = self.assign_h_index_matrix()
    self.n_h = self.calc_nh()

  def create_Ih_table(self, refl_table, weights):
    '''create an Ih_table from the reflection table'''
    #check necessary columns exists in input reflection table
    columns = ['asu_miller_index', 'intensity', 'inverse_scale_factor', 'Esq']
    for col in columns:
      if not col in refl_table.keys():
        assert 0, """Attempting to create an Ih_table object from a reflection
        table with no %s column""" % col
    if len(refl_table) != len(weights):
      assert 0, """Attempting to create an Ih_table object from a reflection
      table and weights list of unequal length."""
    Ih_table = flex.reflection_table()
    for col in columns:
      Ih_table[col] = refl_table[col]
    Ih_table['Ih_values'] = flex.double([0.0] * len(refl_table))
    Ih_table['weights'] = weights
    return Ih_table.select(Ih_table['weights'] != 0.0)

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

  def assign_h_index(self):
    '''assign an index to the sorted reflection table that
       labels each group of unique miller indices'''
    h_index_counter_array = []
    h_index = 0
    h_index_counter = 1
    for i in range(1, len(self.Ih_table)):
      if (self.Ih_table['asu_miller_index'][i] ==
          self.Ih_table['asu_miller_index'][i-1]):
        h_index_counter += 1
      else:
        h_index += 1
        h_index_counter_array.append(h_index_counter)
        h_index_counter = 1
    h_index_counter_array.append(h_index_counter)
    '''calculate the cumulative sum after each h_index group'''
    hsum = 0
    h_index_cumulative_array = [0]
    for n in h_index_counter_array:
      hsum += n
      h_index_cumulative_array.append(hsum)
    return flex.int(h_index_counter_array), flex.int(h_index_cumulative_array)

  def assign_h_index_matrix(self):
    n1 = len(self.Ih_table['asu_miller_index'])
    h_index_mat = sparse.matrix(n1, len(self.h_index_counter_array))
    for i in range(len(self.h_index_cumulative_array)-1):
      col = sparse.matrix_column(n1)
      start_idx = self.h_index_cumulative_array[i]
      for j in range(self.h_index_counter_array[i]):
        col[start_idx+j] = 1
      h_index_mat[:, i] = col
    return h_index_mat

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
    scale_factors = self.Ih_table['inverse_scale_factor']
    gsq = (((scale_factors)**2) * self.Ih_table['weights'])
    sumgsq = gsq * self.h_index_mat
    gI = ((scale_factors * self.Ih_table['intensity']) * self.Ih_table['weights'])
    sumgI = gI * self.h_index_mat
    Ih = sumgI/sumgsq
    self.Ih_table['Ih_values'] = flex.double(
      np.repeat(Ih, self.h_index_counter_array))


class joined_Ih_table(base_Ih_table):
  def __init__(self, datamanagers):
    self.Ih_tables = []
    self.experiments = []
    for dm in datamanagers:
      self.Ih_tables.append(dm.Ih_table)
      self.experiments.append(dm.experiments)
    self.h_idx_count_list = []
    self.h_idx_cumulative_list = []
    self.Ih_table, self.unique_indices = self.determine_all_unique_indices()
    self.assign_h_index_arrays()
    self.h_index_mat = self.assign_h_index_matrix()
    self.h_idx_expand_list = self.assign_h_expand_matrices()
    #finish construction of Ih_table
    self.construct_Ih_table()
    self.n_h = self.calc_nh()

  def determine_all_unique_indices(self):
    '''this function finds the unique reflections across all datasets and
    also writes the sorted asu miller indices to an Ih_table.'''
    s_g_1 = self.experiments[0].crystal.get_space_group()
    for experiment in self.experiments:
      assert experiment.crystal.get_space_group() == s_g_1
    crystal_symmetry = crystal.symmetry(space_group=s_g_1)
    all_miller_indices = []
    for Ih_tab in self.Ih_tables:
      all_miller_indices.extend(list(Ih_tab.Ih_table['asu_miller_index']))

    all_unique_indices = flex.miller_index(list(set(all_miller_indices)))
    unique_miller_set = miller.set(crystal_symmetry=crystal_symmetry,
                                   indices=all_unique_indices)
    unique_permuted = unique_miller_set.sort_permutation(by_value='packed_indices')
    unique_indices = all_unique_indices.select(unique_permuted)

    all_miller_indices = flex.miller_index(all_miller_indices)
    full_miller_set = miller.set(crystal_symmetry=crystal_symmetry,
                                 indices=all_miller_indices)
    all_permuted = full_miller_set.sort_permutation(by_value='packed_indices')
    Ih_table = flex.reflection_table()
    Ih_table['asu_miller_index'] = all_miller_indices.select(all_permuted)
    return Ih_table, unique_indices


  def assign_h_index_arrays(self):
    '''this function determines new h_index counter and cumulative arrays for
    the individual datasets to enable later creation of the h_expander matrices.
    The counter and cumulative arrays for the joint dataset are also created.'''
    for Ih_table in self.Ih_tables:
      miller_idx = Ih_table.Ih_table['asu_miller_index']
      h_idx_count = flex.int([])
      #note: different to single case as need to count the zero instances as well
      for unique_index in self.unique_indices:
        n = (miller_idx == unique_index).count(True)
        h_idx_count.append(n)
      hsum = 0
      h_index_cumulative_array = flex.int([0])
      for n in h_idx_count:
        hsum += n
        h_index_cumulative_array.append(hsum)
      self.h_idx_count_list.append(h_idx_count)
      self.h_idx_cumulative_list.append(h_index_cumulative_array)

    self.h_index_counter_array = flex.int([0]*len(self.h_idx_count_list[0]))
    self.h_index_cumulative_array = flex.int([0]*len(self.h_idx_cumulative_list[0]))
    for h_idx_count, h_index_cumul in zip(self.h_idx_count_list, self.h_idx_cumulative_list):
      self.h_index_counter_array += h_idx_count
      self.h_index_cumulative_array += h_index_cumul

  def assign_h_expand_matrices(self):
    '''this function creates a h_expand matrix for each dataset, so that
    quantites from the individual datasets can easily be combined during
    scaling whilst maintaining the ordered grouping into symmetry
    equivalent reflections.'''
    n_total_refl = self.h_index_cumulative_array[-1]
    h_idx_expand_list = []
    for m, Ih_table in enumerate(self.Ih_tables):
      n_refl = len(Ih_table.Ih_table)
      h_expand_mat = sparse.matrix(n_refl, n_total_refl)
      #delete certain elements to make the idx calculation easy in the loop:
      #shift arrays by -1 for n < m. This will probably end badly.
      for n in range(m):
        del self.h_idx_cumulative_list[n][0] #this element was zero
      last_elements = []
      for n in range(m, len(self.Ih_tables)):
        last_elements.append(self.h_idx_cumulative_list[n][-1])
        del self.h_idx_cumulative_list[n][-1] #remove last element of the rest
      summed_cumulative_arrays = flex.int([0]*len(self.h_idx_cumulative_list[0]))
      for n in range(len(self.Ih_tables)):
        summed_cumulative_arrays += self.h_idx_cumulative_list[n]

      counter = 0
      for i, val in enumerate(self.h_idx_count_list[m]):
        for j in range(val):
          idx = j + summed_cumulative_arrays[i]
          h_expand_mat[counter, idx] = 1
          counter += 1
      h_idx_expand_list.append(h_expand_mat)

      #now put the elements back into the cumulative arrays for the next calc!
      for n in range(m):
        self.h_idx_cumulative_list[n].insert(0, 0)
      counter = 0
      for n in range(m, len(self.Ih_tables)):
        self.h_idx_cumulative_list[n].append(last_elements[counter])
        counter += 1
    return h_idx_expand_list

  def construct_Ih_table(self):
    '''method to allow finishing the construction of the Ih_table now that
    the expand matrices have been determined.'''
    n = self.h_index_cumulative_array[-1]
    intensities = flex.double([0.0]*n)
    scales = flex.double([0.0]*n)
    scaleweights = flex.double([0.0]*n)

    for i, Ih_table in enumerate(self.Ih_tables):
      intensities += Ih_table.Ih_table['intensity'] * self.h_idx_expand_list[i]
      scales += Ih_table.Ih_table['inverse_scale_factor'] * self.h_idx_expand_list[i]
      scaleweights += Ih_table.Ih_table['weights'] * self.h_idx_expand_list[i]

    gsq = (((scales)**2) * scaleweights)
    sumgsq = gsq * self.h_index_mat
    gI = ((scales * intensities) * scaleweights)
    sumgI = gI * self.h_index_mat
    Ih = sumgI/sumgsq
    self.Ih_table['intensity'] = intensities
    self.Ih_table['inverse_scale_factor'] = scales
    self.Ih_table['weights'] = scaleweights
    self.Ih_table['Ih_values'] = flex.double(np.repeat(Ih, self.h_index_counter_array))

  def calc_Ih(self):
    '''calculate the current best estimate for I for each reflection group'''
    scales = flex.double([0.0] * self.h_index_cumulative_array[-1])
    for i, Ih_table in enumerate(self.Ih_tables):
      scales += Ih_table.Ih_table['inverse_scale_factor'] * self.h_idx_expand_list[i]
    gsq = (((scales)**2) * self.Ih_table['weights'])
    sumgsq = gsq * self.h_index_mat
    gI = ((scales * self.Ih_table['intensity']) * self.Ih_table['weights'])
    sumgI = gI * self.h_index_mat
    Ih = sumgI/sumgsq
    self.Ih_table['inverse_scale_factor'] = scales
    self.Ih_table['Ih_values'] = flex.double(np.repeat(Ih, self.h_index_counter_array))
