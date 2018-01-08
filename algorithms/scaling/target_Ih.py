'''
Define the data structures required for a scaling algorithm. This consists of
two distinct classes - a SingleIhTable and JointIhTable, which inherit
properties from IhTableBase. Both SingleIhTable and JointIhTable have a similar
interface, consisting of the public attributes size, h_index_counter_array,
h_index_cumulative_array, h_index_matrix and n_h (the number of memebers
in each group of equivalent reflections), as well as access to the data via
the attributes weights, intensities, inverse_scale_factors, asu_miller_index
and Ih_values. The JointIhTable also has a h_index_expand_list attribute.
'''
import abc
import numpy as np
from dials.array_family import flex
from cctbx import miller, crystal
from scitbx import sparse

class IhTableBase(object):
  '''Base class of methods needed to create the datastructure
  for a scaling algorithm'''
  __metaclass__ = abc.ABCMeta

  def __init__(self, data):
    self._h_index_counter_array = None
    self._h_index_cumulative_array = None
    self._h_index_matrix = None
    self._n_h = None
    self._Ih_table = self._create_Ih_table(data)

  @abc.abstractmethod
  def _create_Ih_table(self, data):
    """Create an Ih_table using the constructor in a subclass"""
    pass

  @property
  def size(self):
    return self._Ih_table.size()

  @property
  def weights(self):
    return self._Ih_table['weights']

  @weights.setter
  def weights(self, new_weights):
    '''method to set new weights'''
    if len(new_weights) != len(self.weights):
      assert 0, """attempting to set a new set of weights of different
      length than previous assignment: was %s, attempting %s""" % (
        len(self.weights), len(new_weights))
    self._Ih_table['weights'] = new_weights

  @property
  def intensities(self):
    return self._Ih_table['intensity']

  @property
  def inverse_scale_factors(self):
    return self._Ih_table['inverse_scale_factor']

  @inverse_scale_factors.setter
  def inverse_scale_factors(self, new_scales):
    '''method to set new scale factors'''
    if len(new_scales) != len(self.inverse_scale_factors):
      assert 0, """attempting to set a new set of scale factors of different
      length than previous assignment: was %s, attempting %s""" % (
        len(self.inverse_scale_factors), len(new_scales))
    self._Ih_table['inverse_scale_factor'] = new_scales

  @property
  def Ih_values(self):
    return self._Ih_table['Ih_values']

  @property
  def asu_miller_index(self):
    return self._Ih_table['asu_miller_index']

  @property
  def h_index_counter_array(self):
    return self._h_index_counter_array

  @property
  def h_index_cumulative_array(self):
    return self._h_index_cumulative_array

  @property
  def h_index_matrix(self):
    return self._h_index_matrix

  @property
  def n_h(self):
    return self._n_h

  def update_aimless_error_model(self, error_params):
    '''method to update the scaling weights using an aimless error model'''
    #note - does this mean we should keep the initial variances?
    #ok for now if you don't update the error model more than once?
    sigmaprime = (((1.0/self.weights) + ((error_params[1] * self.intensities)**2)
                  )**0.5) * error_params[0]
    self.weights = 1.0/(sigmaprime**2)

  def select(self, selection):
    ''''selects a subset of the data and recalculates h_index_arrays/matrix,
    before returning self (to act like a flex selection operation).'''
    self._Ih_table = self._Ih_table.select(selection)
    (self._h_index_counter_array, self._h_index_cumulative_array
    ) = self._assign_h_index(self.asu_miller_index)
    self._h_index_matrix = self._assign_h_index_matrix(
      self.h_index_counter_array, self.h_index_cumulative_array)
    self._n_h = self._calc_nh(self.h_index_counter_array)
    return self

  @staticmethod
  def _assign_h_index(asu_miller_index):
    '''assign an index to the Ih table that
       labels each group of unique miller indices'''
    h_index_counter_array = []
    h_index = 0
    h_index_counter = 1
    for i in range(1, len(asu_miller_index)):
      if asu_miller_index[i] == asu_miller_index[i-1]:
        h_index_counter += 1
      else:
        h_index += 1
        h_index_counter_array.append(h_index_counter)
        h_index_counter = 1
    h_index_counter_array.append(h_index_counter)
    #now calculate the cumulative sum after each h_index group
    hsum = 0
    h_index_cumulative_array = [0]
    for n in h_index_counter_array:
      hsum += n
      h_index_cumulative_array.append(hsum)
    return flex.int(h_index_counter_array), flex.int(h_index_cumulative_array)

  @staticmethod
  def _assign_h_index_matrix(h_idx_count_arr, h_idx_cumul_arr):
    '''assign a h_index matrix to allow fast summation over groups of
    symmetry equivalent reflections in the Ih_table'''
    n1 = h_idx_cumul_arr[-1]
    h_index_matrix = sparse.matrix(n1, len(h_idx_count_arr))
    for i in range(len(h_idx_cumul_arr)-1):
      col = sparse.matrix_column(n1)
      start_idx = h_idx_cumul_arr[i]
      for j in range(h_idx_count_arr[i]):
        col[start_idx+j] = 1
      h_index_matrix[:, i] = col
    return h_index_matrix

  @staticmethod
  def _calc_nh(h_index_counter_array):
    '''returns a vector of len(reflections) with the number of members of
    each h group'''
    n_h = flex.double([])
    for i in h_index_counter_array:
      n_h.extend(flex.double([i]*i))
    return n_h


class SingleIhTable(IhTableBase):
  '''Class to create an Ih_table. This is the default
  data structure used for scaling a single sweep.'''
  def __init__(self, reflection_table, weighting):
    super(SingleIhTable, self).__init__([reflection_table, weighting])
    (self._h_index_counter_array, self._h_index_cumulative_array
    ) = self._assign_h_index(self.asu_miller_index)
    self._h_index_matrix = self._assign_h_index_matrix(
      self.h_index_counter_array, self._h_index_cumulative_array)
    self._n_h = self._calc_nh(self.h_index_counter_array)
    self.calc_Ih() #calculate a first estimate of Ih

  def _create_Ih_table(self, data):
    '''create an Ih_table from the reflection table'''
    (refl_table, weights) = data
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

  def calc_Ih(self):
    '''calculate the current best estimate for I for each reflection group'''
    scale_factors = self.inverse_scale_factors
    gsq = (((scale_factors)**2) * self.weights)
    sumgsq = gsq * self.h_index_matrix
    gI = ((scale_factors * self.intensities) * self.weights)
    sumgI = gI * self.h_index_matrix
    Ih = sumgI/sumgsq
    self._Ih_table['Ih_values'] = flex.double(np.repeat(Ih, self.h_index_counter_array))


class JointIhTable(IhTableBase):
  '''Class to expand the datastructure for scaling multiple
  datasets together.'''
  def __init__(self, datamanagers):
    self._h_index_expand_list = None
    super(JointIhTable, self).__init__(data=datamanagers)
    self._n_h = self._calc_nh(self.h_index_counter_array)

  @property
  def h_index_expand_list(self):
    return self._h_index_expand_list

  def _create_Ih_table(self, data):
    '''construct a single Ih_table for the combined reflections, using the
    Ih_tables from the individual datamanagers.'''
    self._Ih_tables = []
    self._experiments = []
    for datamanager in data:
      self._Ih_tables.append(datamanager.Ih_table)
      self._experiments.append(datamanager.experiments)
    self._h_idx_count_list = []
    self._h_idx_cumulative_list = []
    Ih_table, self._unique_indices = self._determine_all_unique_indices()
    self._assign_h_index_arrays()
    self._h_index_matrix = self._assign_h_index_matrix(
      self.h_index_counter_array, self.h_index_cumulative_array)
    self._h_index_expand_list = self._assign_h_expand_matrices()
    #finish construction of Ih_table
    Ih_table = self._complete_Ih_table(Ih_table)
    return Ih_table

  def calc_Ih(self):
    '''calculate the current best estimate for I for each reflection group'''
    scales = flex.double([0.0] * self.h_index_cumulative_array[-1])
    for i, Ih_table in enumerate(self._Ih_tables):
      scales += Ih_table.inverse_scale_factors * self.h_index_expand_list[i]
    gsq = (((scales)**2) * self.weights)
    sumgsq = gsq * self.h_index_matrix
    gI = ((scales * self.intensities) * self.weights)
    sumgI = gI * self.h_index_matrix
    Ih = sumgI/sumgsq
    self.inverse_scale_factors = scales
    self._Ih_table['Ih_values'] = flex.double(np.repeat(Ih, self.h_index_counter_array))

  def _determine_all_unique_indices(self):
    '''this function finds the unique reflections across all datasets and
    also writes the sorted asu miller indices to an Ih_table.'''
    s_g_1 = self._experiments[0].crystal.get_space_group()
    for experiment in self._experiments:
      assert experiment.crystal.get_space_group() == s_g_1
    crystal_symmetry = crystal.symmetry(space_group=s_g_1)
    all_miller_indices = []
    for Ih_tab in self._Ih_tables:
      all_miller_indices.extend(list(Ih_tab.asu_miller_index))
    #find the ordered set of unique indices
    all_unique_indices = flex.miller_index(list(set(all_miller_indices)))
    unique_miller_set = miller.set(crystal_symmetry=crystal_symmetry,
                                   indices=all_unique_indices)
    unique_permuted = unique_miller_set.sort_permutation(by_value='packed_indices')
    unique_indices = all_unique_indices.select(unique_permuted)
    #find the ordered set of all indices
    all_miller_indices = flex.miller_index(all_miller_indices)
    full_miller_set = miller.set(crystal_symmetry=crystal_symmetry,
                                 indices=all_miller_indices)
    all_permuted = full_miller_set.sort_permutation(by_value='packed_indices')
    Ih_table = flex.reflection_table()
    Ih_table['asu_miller_index'] = all_miller_indices.select(all_permuted)
    return Ih_table, unique_indices

  def _assign_h_index_arrays(self):
    '''this function determines new h_index counter and cumulative arrays for
    the individual datasets to enable later creation of the h_expander matrices.
    The counter and cumulative arrays for the joint dataset are also created.'''
    for Ih_table in self._Ih_tables:
      miller_idx = Ih_table.asu_miller_index
      h_idx_count = flex.int([])
      #note: different to single case as need to count the zero instances as well
      for unique_index in self._unique_indices:
        n = (miller_idx == unique_index).count(True)
        h_idx_count.append(n)
      hsum = 0
      h_index_cumulative_array = flex.int([0])
      for n in h_idx_count:
        hsum += n
        h_index_cumulative_array.append(hsum)
      self._h_idx_count_list.append(h_idx_count)
      self._h_idx_cumulative_list.append(h_index_cumulative_array)
    #now calculate the cumulative/counter arrays for the joint dataset.
    self._h_index_counter_array = flex.int([0]*len(self._h_idx_count_list[0]))
    self._h_index_cumulative_array = flex.int([0]*len(self._h_idx_cumulative_list[0]))
    for h_idx_count, h_index_cumul in zip(self._h_idx_count_list, self._h_idx_cumulative_list):
      self._h_index_counter_array += h_idx_count
      self._h_index_cumulative_array += h_index_cumul

  def _assign_h_expand_matrices(self):
    '''this function creates a h_expand matrix for each dataset, so that
    quantites from the individual datasets can easily be combined during
    scaling whilst maintaining the ordered grouping into symmetry
    equivalent reflections.'''
    n_total_refl = self.h_index_cumulative_array[-1]
    h_idx_expand_list = []
    for m, Ih_table in enumerate(self._Ih_tables):
      n_refl = Ih_table.size
      h_expand_mat = sparse.matrix(n_refl, n_total_refl)
      #delete certain elements to make the idx calculation easy in the loop:
      #shift arrays by -1 for n < m.
      for n in range(m):
        del self._h_idx_cumulative_list[n][0] #this element was zero
      last_elements = []
      for n in range(m, len(self._Ih_tables)):
        last_elements.append(self._h_idx_cumulative_list[n][-1])
        del self._h_idx_cumulative_list[n][-1] #remove last element of the rest
      summed_cumulative_arrays = flex.int([0]*len(self._h_idx_cumulative_list[0]))
      for n in range(len(self._Ih_tables)):
        summed_cumulative_arrays += self._h_idx_cumulative_list[n]
      #now easily construct the h_expand_matrix
      counter = 0
      for i, val in enumerate(self._h_idx_count_list[m]):
        for j in range(val):
          idx = j + summed_cumulative_arrays[i]
          h_expand_mat[counter, idx] = 1
          counter += 1
      h_idx_expand_list.append(h_expand_mat)
      #now put the elements back into the cumulative arrays for the next calc!
      for n in range(m):
        self._h_idx_cumulative_list[n].insert(0, 0)
      counter = 0
      for n in range(m, len(self._Ih_tables)):
        self._h_idx_cumulative_list[n].append(last_elements[counter])
        counter += 1
    return h_idx_expand_list

  def _complete_Ih_table(self, Ih_table):
    '''method to allow finishing the construction of the Ih_table now that
    the expand matrices have been determined.'''
    n = self.h_index_cumulative_array[-1]
    intensities = flex.double([0.0]*n)
    scales = flex.double([0.0]*n)
    scaleweights = flex.double([0.0]*n)
    for i, Ih_tab in enumerate(self._Ih_tables):
      intensities += Ih_tab.intensities * self.h_index_expand_list[i]
      scales += Ih_tab.inverse_scale_factors * self.h_index_expand_list[i]
      scaleweights += Ih_tab.weights * self.h_index_expand_list[i]
    gsq = (((scales)**2) * scaleweights)
    sumgsq = gsq * self.h_index_matrix
    gI = ((scales * intensities) * scaleweights)
    sumgI = gI * self.h_index_matrix
    Ih = sumgI/sumgsq
    Ih_table['intensity'] = intensities
    Ih_table['inverse_scale_factor'] = scales
    Ih_table['weights'] = scaleweights
    Ih_table['Ih_values'] = flex.double(np.repeat(Ih, self.h_index_counter_array))
    return Ih_table
