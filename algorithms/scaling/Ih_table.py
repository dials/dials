"""
Classes that define the sorted datastructures needed for scaling.

The basic design of the data structure, called the Ih_table, is
similar to a reflection_table but is initialised with a sorted
reflection table, where reflections are grouped by their asu
miller index.
In addition to the table, a set of sparse matrices are defined;
the h_index_matrix, which is used to efficiently calculate sums
over groups of reflections, and the h_expand_matrix, which is
used to in a similar way to numpy.repeat, to repeat summation
properties so that calculations are fully vectorised.
Data structures are defined for scaling single and multiple
datasets simulatenously - the JointIhTable has further matrices
to keep track of the locations of equivalent reflections across
datasets and still perform fully vectorised calculations.
Access to the data is given through the attributes weights,
intensities, inverse_scale_factors, asu_miller_index and Ih_values.
"""
import abc
import numpy as np
from dials.array_family import flex
from cctbx import miller, crystal
from scitbx import sparse

class IhTableBase(object):
  """Base class that defines the interface of the datastructure."""

  __metaclass__ = abc.ABCMeta

  def __init__(self, data):
    self._h_index_matrix = None
    self._h_expand_matrix = None
    self._Ih_table = self._create_Ih_table(data)
    self._n_h = None

  @abc.abstractmethod
  def _create_Ih_table(self, data):
    """Create an Ih_table using the constructor in a subclass."""
    pass

  @abc.abstractmethod
  def _assign_h_matrices(self, asu_miller_index):
    """Assign the h_index and h_expand matrices."""
    pass

  @abc.abstractmethod
  def calc_Ih(self):
    """Calculate the current best estimate for I for each reflection group."""
    pass

  @property
  def size(self):
    """Return the length of the stored Ih_table (a reflection table)."""
    return self._Ih_table.size()

  @property
  def weights(self):
    """The weights that will be used in scaling."""
    return self._Ih_table['weights']

  @weights.setter
  def weights(self, new_weights):
    if new_weights.size() != self.weights.size():
      assert 0, '''attempting to set a new set of weights of different
      length than previous assignment: was %s, attempting %s''' % (
        self.weights.size(), new_weights.size())
    self._Ih_table['weights'] = new_weights

  @property
  def variances(self):
    """The initial variances of the reflections."""
    return self._Ih_table['variance']

  @property
  def intensities(self):
    """The unscaled reflection intensities."""
    return self._Ih_table['intensity']

  @property
  def inverse_scale_factors(self):
    """"The inverse scale factors of the reflections."""
    return self._Ih_table['inverse_scale_factor']

  @inverse_scale_factors.setter
  def inverse_scale_factors(self, new_scales):
    if new_scales.size() != self.inverse_scale_factors.size():
      assert 0, """attempting to set a new set of scale factors of different
      length than previous assignment: was %s, attempting %s""" % (
        self.inverse_scale_factors.size(), new_scales.size())
    self._Ih_table['inverse_scale_factor'] = new_scales

  @property
  def Ih_values(self):
    """The estimated intensities of symmetry equivalent reflections."""
    return self._Ih_table['Ih_values']

  @property
  def asu_miller_index(self):
    """The asymmetric miller indices."""
    return self._Ih_table['asu_miller_index']

  @property
  def miller_index(self):
    """The miller indices."""
    return self._Ih_table['miller_index']

  @property
  def h_index_matrix(self):
    """A sparse matrix to perform sums over groups of unique reflections.

    Given a flex array a, the sum over unique reflections is given
    by a * h_index_matrix. h_index_matrix is an n_reflections x
    n_symmetry_unique_groups matrix. The only nonzero elements in the nth
    column have values of 1.0, in rows corresponding to the positions of the
    nth symmetry equivalent group in the flex array a."""
    return self._h_index_matrix

  @property
  def h_expand_matrix(self):
    """A sparse matrix to expand out a property obtained by a sum over
    unique reflections.

    For example, wgI_sum = sum_h weights * scales * intensity, is a vector
    of length n_symmetry_unique_groups. wgI * h_expand_matrix is then
    a vector of length n_reflections, containing the wgI_sum corresponding
    to the symmetry group of each reflection. This can then be used for
    vectorised calculations. h_expand_matrix is the transpose of the
    h_index_matrix."""
    return self._h_expand_matrix

  @property
  def n_h(self):
    """A vector of length n_refl, containing the number of reflections in
    the respective symmetry group.
    
    Not calculated by default, as only needed for certain calculations."""
    return self._n_h

  def calc_nh(self):
    """Calculate the n_h vector."""
    return ((flex.double([1.0] * self.size) * self.h_index_matrix)
      * self.h_expand_matrix)

  def update_error_model(self, error_params):
    """Update the scaling weights based on an error model."""
    sigmaprime = (((self.variances) + ((error_params[1] * self.intensities)**2)
                  )**0.5) * error_params[0]
    self.weights = 1.0/(sigmaprime**2)

  def apply_tukey_biweighting(self):
    """Apply a tukey biweighting scheme for the scaling weights."""
    z_score = flex.double([])
    zmax = 6.0
    for i, _ in enumerate(self.h_index_counter_array):
      h_idx_cumul = self.h_index_cumulative_array[i:i+2]
      Ihls = self.intensities[h_idx_cumul[0]:h_idx_cumul[1]]
      var = self.variances[h_idx_cumul[0]:h_idx_cumul[1]]
      med = np.median(Ihls)
      sigma = max([np.median(var**0.5), np.median(Ihls - med)])
      z = (Ihls - med) / sigma
      z_score.extend(z)
    tukey_weights = (1.0 - ((z_score/zmax)**2))**2
    sel = tukey_weights < 0.0
    tukey_weights.set_selected(sel, 0.0)
    self.weights = tukey_weights

  def select(self, selection):
    """Select a subset of the data and recalculate h_index_matrices,
    before returning self (to act like a flex selection operation)."""
    self._Ih_table = self._Ih_table.select(selection)
    self._h_index_matrix, self._h_expand_matrix = self._assign_h_matrices(
      self.asu_miller_index)
    return self

class SingleIhTable(IhTableBase):
  """Class to create an Ih_table. This is the default
  data structure used for scaling a single sweep."""

  def __init__(self, reflection_table, weights=None):
    self._nonzero_weights = None
    super(SingleIhTable, self).__init__([reflection_table, weights])
    self._h_index_matrix, self._h_expand_matrix = self._assign_h_matrices(
      self.asu_miller_index)
    if not 'Ih_values' in reflection_table.keys():
      self.calc_Ih() #calculate a first estimate of Ih

  @property
  def nonzero_weights(self):
    """Retain selection array relative to input reflection table, to use
    for referring outliers back to initial input."""
    return self._nonzero_weights

  def _create_Ih_table(self, data):
    """Create an Ih_table from the reflection table and optionally weights.

    It is assumed that the refl_table is already sorted by asu_miller_index."""
    (refl_table, weights) = data
    Ih_table = flex.reflection_table()
    for col in ['asu_miller_index', 'intensity', 'inverse_scale_factor',
      'variance', 'miller_index']:
      if not col in refl_table.keys():
        assert 0, """Attempting to create an Ih_table object from a reflection
        table with no %s column""" % col
      Ih_table[col] = refl_table[col]
    if 'Ih_values' in refl_table.keys():
      Ih_table['Ih_values'] = refl_table['Ih_values']
    else:
      Ih_table['Ih_values'] = flex.double([0.0] * refl_table.size())
    if weights:
      if refl_table.size() != weights.size():
        assert 0, """Attempting to create an Ih_table object from a reflection
        table and weights list of unequal length."""
      Ih_table['weights'] = weights
      nonzero_weights_sel = Ih_table['weights'] != 0.0
      Ih_table = Ih_table.select(nonzero_weights_sel)
    else:
      nonzero_weights_sel = ~(
        refl_table.get_flags(refl_table.flags.user_excluded_in_scaling)
        | refl_table.get_flags(refl_table.flags.excluded_for_scaling))
      Ih_table = Ih_table.select(nonzero_weights_sel)
      Ih_table['weights'] = 1.0/Ih_table['variance']
    self._nonzero_weights = nonzero_weights_sel
    return Ih_table

  def _assign_h_matrices(self, asu_miller_index):
    """Assign the h_index and h_expand matrices."""
    n_refl = asu_miller_index.size()
    n_unique_groups = len(set(asu_miller_index))
    h_index_matrix = sparse.matrix(n_refl, n_unique_groups)
    group_idx = 0
    previous = asu_miller_index[0]
    for refl_idx, asu_idx in enumerate(asu_miller_index):
      if asu_idx != previous:
        group_idx += 1
      h_index_matrix[refl_idx, group_idx] = 1.0
      previous = asu_idx
    h_expand_matrix = h_index_matrix.transpose()
    return h_index_matrix, h_expand_matrix

  def calc_Ih(self):
    """Calculate the current best estimate for I for each reflection group."""
    scale_factors = self.inverse_scale_factors
    gsq = (((scale_factors)**2) * self.weights)
    sumgsq = gsq * self.h_index_matrix
    gI = ((scale_factors * self.intensities) * self.weights)
    sumgI = gI * self.h_index_matrix
    Ih = sumgI/sumgsq
    self._Ih_table['Ih_values'] = Ih * self.h_expand_matrix

class JointIhTable(IhTableBase):
  """Class to expand the datastructure for scaling multiple
  datasets together."""

  def __init__(self, single_scalers):
    self._Ih_tables = []
    self._experiments = []
    self._unique_indices = None
    super(JointIhTable, self).__init__(data=single_scalers)

  def _create_Ih_table(self, data):
    """Construct a single Ih_table for the combined reflections, using the
    Ih_tables from the individual scalers.

    The data is stored differently to a single Ih_Table - the data here is
    not sorted by miller index, but is formed by extending each of the
    individual sorted tables. This does not matter, as the order is encoded
    within the h_index, h_expand matrices. Keeping the data in extended order
    allows a quicker calculation of Ih."""
    for scaler in data:
      self._Ih_tables.append(scaler.Ih_table)
      self._experiments.append(scaler.experiments)
    self._unique_indices, all_indices = (
      self._determine_all_unique_indices())
    self._h_index_matrix, self._h_expand_matrix = self._assign_h_matrices(
      all_indices)
    # Now join together the datasets to make the Ih table
    Ih_table = flex.reflection_table()
    for Ih_tab in self._Ih_tables:
      Ih_table.extend(Ih_tab._Ih_table)
    scales = Ih_table['inverse_scale_factor']
    scaleweights = Ih_table['weights']
    intensities = Ih_table['intensity']
    sumgsq = (((scales)**2) * scaleweights) * self.h_index_matrix
    sumgI = ((scales * intensities) * scaleweights) * self.h_index_matrix
    Ih_table['Ih_values'] = (sumgI/sumgsq) * self.h_expand_matrix
    return Ih_table

  def _determine_all_unique_indices(self):
    """Determine the ordered set of unique indices and the ordered set
    of all asu miller indices across all datasets."""
    s_g_1 = self._experiments[0].crystal.get_space_group()
    for experiment in self._experiments:
      assert experiment.crystal.get_space_group() == s_g_1
    crystal_symmetry = crystal.symmetry(space_group=s_g_1)
    all_miller_indices = []
    for Ih_tab in self._Ih_tables:
      all_miller_indices.extend(list(Ih_tab.asu_miller_index))
    # Find the ordered set of unique indices
    all_unique_indices = flex.miller_index(list(set(all_miller_indices)))
    unique_ms = miller.set(crystal_symmetry=crystal_symmetry,
      indices=all_unique_indices)
    unique_permuted = unique_ms.sort_permutation(by_value='packed_indices')
    unique_indices = all_unique_indices.select(unique_permuted)
    # Find the ordered set of all indices
    all_miller_indices = flex.miller_index(all_miller_indices)
    full_miller_set = miller.set(crystal_symmetry=crystal_symmetry,
      indices=all_miller_indices)
    all_permuted = full_miller_set.sort_permutation(by_value='packed_indices')
    all_indices = all_miller_indices.select(all_permuted)
    return unique_indices, all_indices

  def _assign_h_matrices(self, asu_miller_index):
    """Assign the h_index and h_expand matrices."""
    n_unique_groups = self._unique_indices.size()
    total_h_idx_matrix = sparse.matrix(asu_miller_index.size(), n_unique_groups)
    h_idx_matrix_row = 0
    crystal_symmetry = crystal.symmetry(
      space_group=self._experiments[0].crystal.get_space_group())
    for Ih_table in self._Ih_tables:
      h_expand_matrix = sparse.matrix(Ih_table.h_index_matrix.n_cols,
        n_unique_groups)
      unique_indices = flex.miller_index(list(set(Ih_table.asu_miller_index)))
      unique_miller_set = miller.set(crystal_symmetry=crystal_symmetry,
        indices=unique_indices)
      unique_indices_sorted = unique_indices.select(
        unique_miller_set.sort_permutation(by_value='packed_indices'))
      # Extend by one to allow a quick loop below without bounds check.
      unique_indices_sorted.extend(flex.miller_index([(0, 0, 0)]))
      indiv_group_idx = 0
      for overall_group_idx, asu_idx in enumerate(self._unique_indices):
        if asu_idx == unique_indices_sorted[indiv_group_idx]:
          h_expand_matrix[indiv_group_idx, overall_group_idx] = 1.0
          indiv_group_idx += 1
      h_index_component = Ih_table.h_index_matrix * h_expand_matrix
      total_h_idx_matrix.assign_block(h_index_component, h_idx_matrix_row, 0)
      h_idx_matrix_row += Ih_table.h_index_matrix.n_rows
    total_h_expand_matrix = total_h_idx_matrix.transpose()
    return total_h_idx_matrix, total_h_expand_matrix

  def calc_Ih(self):
    """Calculate the current best estimate for I for each reflection group."""
    scales = flex.double([])
    for Ih_table in self._Ih_tables:
      scales.extend(Ih_table.inverse_scale_factors)
    sumgsq = (((scales)**2) * self.weights) * self.h_index_matrix
    sumgI = ((scales * self.intensities) * self.weights) * self.h_index_matrix
    Ih = sumgI/sumgsq
    self.inverse_scale_factors = scales
    self._Ih_table['Ih_values'] = Ih * self.h_expand_matrix

  def update_weights_from_error_models(self):
    """Reset weights, to be used after the weights of individual Ih tables
    have been modified."""
    weights = flex.double()
    for Ih_tab in self._Ih_tables:
      weights.extend(Ih_tab.weights)
    self.weights = weights
