"""
Classes that define the datastructures needed for scaling.

The basic design of the data structure, called the Ih_table, is
similar to a reflection_table. The asu miller index is calculated
and set of sparse matrices are defined; the h_index_matrix, which is used
to efficiently calculate sums over groups of reflections, and the
h_expand_matrix, which is used to in a similar way to numpy.repeat, to
repeat summation properties so that calculations are fully vectorised.
Data structures are defined for scaling single and multiple
datasets simulatenously - the JointIhTable h_index_matrix keeps track of
the locations of equivalent reflections across datasets and still perform
fully vectorised calculations. Access to the data is given through the
attributes weights, intensities, inverse_scale_factors, asu_miller_index
and Ih_values.
"""
import logging
from libtbx.containers import OrderedSet
from cctbx import miller, crystal
from scitbx import sparse
from dials.array_family import flex
from dials.algorithms.scaling.weighting import get_weighting_scheme
logger = logging.getLogger('dials')


class SingleIhTable(object):

  """
  A datastructure for storing a reflection table, with additional matrices for
  efficiently performing sums over groups of equivalent reflections.
  """

  def __init__(self, IhTable, h_index_matrix, space_group, weighting_scheme,
    nonzero_weights_isel):
    self._Ih_table = IhTable
    self.space_group = space_group
    self._h_index_matrix = h_index_matrix
    self._h_expand_matrix = h_index_matrix.transpose()
    self._nonzero_weights = nonzero_weights_isel
    self.derivatives = None
    self._n_h = None
    self.weighting_scheme = get_weighting_scheme(self, weighting_scheme)
    self.weighting_scheme.calculate_initial_weights()
    self.calc_Ih()

  def update_error_model(self, error_params):
    """Update the scaling weights based on an error model."""
    sigmaprime = (((self.variances) + ((error_params[1] * self.intensities)**2)
                  )**0.5) * error_params[0]
    self.weights = 1.0/(sigmaprime**2)

  @property
  def nonzero_weights(self):
    """Retain selection array relative to input reflection table, to use
    for referring outliers back to initial input."""
    return self._nonzero_weights

  def select(self, sel):
    """Select a subset of the """
    self._Ih_table = self._Ih_table.select(sel)
    h_idx_T = self._h_index_matrix.transpose()
    h_idx_sel = h_idx_T.select_columns(sel.iselection())
    reduced_h_idx = h_idx_sel.transpose()
    nz_col_sel = flex.bool(reduced_h_idx.n_cols, True)
    for i, col in enumerate(reduced_h_idx.cols()):
      if col.non_zeroes == 0:
        nz_col_sel[i] = False
    self._h_index_matrix = reduced_h_idx.select_columns(nz_col_sel.iselection())
    self._h_expand_matrix = self._h_index_matrix.transpose()
    self._nonzero_weights = self._nonzero_weights.select(sel)

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
    if new_weights.size() != self.size:
      assert 0, '''attempting to set a new set of weights of different
      length than previous assignment: was %s, attempting %s''' % (
        self.size, new_weights.size())
    self._Ih_table['weights'] = new_weights

  def update_weights(self):
    """Apply an iterative weighting scheme."""
    self.weighting_scheme.apply_iterative_weights()

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
    if new_scales.size() != self.size:
      assert 0, """attempting to set a new set of scale factors of different
      length than previous assignment: was %s, attempting %s""" % (
        self.inverse_scale_factors.size(), new_scales.size())
    else:
      self._Ih_table['inverse_scale_factor'] = new_scales

  @property
  def Ih_values(self):
    """The estimated intensities of symmetry equivalent reflections."""
    return self._Ih_table['Ih_values']

  @property
  def Ih_table(self):
    """A reflection table of all the data stored by the class."""
    return self._Ih_table

  @Ih_table.setter
  def Ih_table(self, Ih_Table):
    """A reflection table of all the data stored by the class."""
    self._Ih_table = Ih_Table

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

  @h_index_matrix.setter
  def h_index_matrix(self, new_matrix):
    assert new_matrix.n_rows == self.size
    self._h_index_matrix = new_matrix

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

  @h_expand_matrix.setter
  def h_expand_matrix(self, new_matrix):
    assert new_matrix.n_cols == self.size
    self._h_expand_matrix = new_matrix

  @property
  def n_h(self):
    """A vector of length n_refl, containing the number of reflections in
    the respective symmetry group.

    Not calculated by default, as only needed for certain calculations."""
    return self._n_h

  def calc_nh(self):
    """Calculate the n_h vector."""
    self._n_h = ((flex.double(self.size, 1.0) * self.h_index_matrix)
      * self.h_expand_matrix)

  def calc_Ih(self):
    """Calculate the current best estimate for I for each reflection group."""
    scale_factors = self.inverse_scale_factors
    gsq = (((scale_factors)**2) * self.weights)
    sumgsq = gsq * self.h_index_matrix
    gI = ((scale_factors * self.intensities) * self.weights)
    sumgI = gI * self.h_index_matrix
    Ih = sumgI/sumgsq
    self._Ih_table['Ih_values'] = Ih * self.h_expand_matrix

  def get_sorted_asu_indices(self):
    """Return the sorted asu indices and the permutation selection."""
    crystal_symmetry = crystal.symmetry(space_group=self.space_group)
    miller_set_in_asu = miller.set(crystal_symmetry=crystal_symmetry,
      indices=self.asu_miller_index, anomalous_flag=False)
    permuted = miller_set_in_asu.sort_permutation(by_value='packed_indices')
    sorted_asu_miller_index = self.asu_miller_index.select(permuted)
    return sorted_asu_miller_index, permuted

  def set_Ih_values_to_target(self, target_Ih_table):
    """Given an Ih table as a target, the common reflections across the tables
    are determined and the Ih_values are set to those of the target. If no
    matching reflection is found, the Ih value is set to zero."""
    target_asu_Ih_dict = dict(zip(target_Ih_table.asu_miller_index,
      target_Ih_table.Ih_values))
    new_Ih_values = flex.double(self.size, 0.0)
    location_in_unscaled_array = 0
    sorted_asu_indices, permuted = self.get_sorted_asu_indices()
    for j, miller_idx in enumerate(OrderedSet(sorted_asu_indices)):
      n_in_group = self.h_index_matrix.col(j).non_zeroes
      if miller_idx in target_asu_Ih_dict:
        i = location_in_unscaled_array
        new_Ih_values.set_selected(flex.size_t(range(i, i + n_in_group)),
          flex.double(n_in_group, target_asu_Ih_dict[miller_idx]))
      location_in_unscaled_array += n_in_group
    self.Ih_values.set_selected(permuted, new_Ih_values)


class IhTable(object):
  """
  A organisational datastructure for storing reflection table data that is
  split into blocks to allow multiprocessing. The data is split by asu miller
  index, i.e. groups of symmetry equivalent reflections are stored in one
  block, even if they originate from different datasets.

  The blocked_data_list attribute is a list of SingleIhTable instances, and
  the blocked_selection_list contains information about the origin of the
  reflections, so that these can be matched back to the input data.

  The preferred initialisation is to pass in a list of (reflection_table,
  selection) tuples - where the reflection_table is the complete dataset (
  including and 'bad' reflections), and selection is a flex.bool array to
  select the data from this table to be used for minimisation. Also required
  is a space group. The number of blocks is controlled by the n_blocks
  parameter, which is one by default, and a weighting scheme can optionally be
  specificed - if none then inverse variance weights are applied."""

  id_ = "IhTable"

  def __init__(self, refl_and_sel_list, space_group, n_blocks=1, weighting_scheme=None):
    self.space_group = space_group
    self.weighting_scheme = weighting_scheme
    self._n_datasets = len(refl_and_sel_list)
    self._blocked_data_list = []
    self._blocked_selection_list = []
    self._free_Ih_table = None
    self._free_set_sel = flex.bool([])
    self._permutation_selection = None
    joint_refl_table, joint_nonzero_weights = self._create_joint_structures(
      refl_and_sel_list)
    joint_refl_table = self._map_indices_to_asu(joint_refl_table)
    self._split_data_into_blocks(joint_refl_table, joint_nonzero_weights,
      n_blocks)
    self.sort_by_dataset_id()
    self._size = joint_refl_table.size()

  def update_error_model(self, error_params):
    """Update the error model in the blocks."""
    for block in self.blocked_data_list:
      block.update_error_model(error_params)

  @property
  def size(self):
    """The total number of reflections in all blocks."""
    return self._size

  @property
  def n_datasets(self):
    """The number of datasets used to construct the Ih_table blocks."""
    return self._n_datasets

  @property
  def blocked_data_list(self):
    """A list of SingleIhTables"""
    return self._blocked_data_list

  @property
  def blocked_selection_list(self):
    """A list of selections relating the Single Ih tables to the initial table."""
    return self._blocked_selection_list

  @property
  def free_Ih_table(self):
    """An Ih table of a subset of the data for a 'free' set in refinement."""
    return self._free_Ih_table

  @property
  def free_set_sel(self):
    """The selection used to split the data into a free and work set."""
    return self._free_set_sel

  @staticmethod
  def _create_joint_structures(refl_and_sel_list):
    """Construct a single Ih_table for the combined reflections, using a list
    of individual Ih_tables.

    The data table of the JointIhTable is formed by extending each of the
    individual tables. The grouping is encoded within the h_index matrix.
    Keeping the data in extended order allows a quicker calculation of Ih."""
    total_Ih_table = flex.reflection_table()
    total_nonzero_weights_sel = flex.size_t()
    for i, (refl_table, sel) in enumerate(refl_and_sel_list):
      Ih_table = flex.reflection_table()
      for col in ['intensity', 'inverse_scale_factor', 'variance', 'miller_index']:
        if not col in refl_table.keys():
          assert 0, """Attempting to create an Ih_table object from a reflection
          table with no %s column""" % col
        Ih_table[col] = refl_table[col]
      if 'Ih_values' in refl_table.keys():
        Ih_table['Ih_values'] = refl_table['Ih_values']
      else:
        Ih_table['Ih_values'] = flex.double(refl_table.size(), 0.0)
      Ih_table['dataset_id'] = flex.int(refl_table.size(), i)
      nonzero_weights_sel = ~(refl_table.get_flags(refl_table.flags.bad_for_scaling,
        all=False))
      if sel:
        nonzero_weights_sel = nonzero_weights_sel & sel
      Ih_table = Ih_table.select(nonzero_weights_sel)
      total_Ih_table.extend(Ih_table)
      total_nonzero_weights_sel.extend(nonzero_weights_sel.iselection())
    return total_Ih_table, total_nonzero_weights_sel

  def _map_indices_to_asu(self, reflection_table):
    """Map the indices to the asymmetric unit."""
    crystal_symmetry = crystal.symmetry(space_group=self.space_group)
    miller_set = miller.set(crystal_symmetry=crystal_symmetry,
      indices=reflection_table['miller_index'], anomalous_flag=False)
    miller_set_in_asu = miller_set.map_to_asu()
    reflection_table['asu_miller_index'] = miller_set_in_asu.indices()
    return reflection_table

  def get_sorted_asu_indices(self, reflection_table):
    """Return the sorted asu indices and the permutation selection."""
    crystal_symmetry = crystal.symmetry(space_group=self.space_group)
    miller_set_in_asu = miller.set(crystal_symmetry=crystal_symmetry,
      indices=reflection_table['asu_miller_index'], anomalous_flag=False)
    permuted = miller_set_in_asu.sort_permutation(by_value='packed_indices')
    sorted_asu_miller_index = reflection_table['asu_miller_index'].select(permuted)
    return sorted_asu_miller_index, permuted

  def _split_data_into_blocks(self, reflection_table, nonzero_weights, n_blocks):
    """Assign the h_index and h_expand matrices."""
    asu_miller_index, self._permutation_selection = \
      self.get_sorted_asu_indices(reflection_table)
    n_refl = asu_miller_index.size()
    n_unique_groups = len(set(asu_miller_index))

    reflection_table = reflection_table.select(self._permutation_selection)
    nonzero_weights = nonzero_weights.select(self._permutation_selection)
    n_blocks = min(n_unique_groups, n_blocks)
    group_boundaries = [int(i*n_unique_groups/n_blocks) for i in range(n_blocks)]
    group_boundaries.append(n_unique_groups)
    total_refl_group_idx = 0 # use to count iteration over all groups
    refl_in_group_idx = 0
    refl_group_idx = 0
    block_idx = 0
    previous = asu_miller_index[0]
    prev_refl_idx = 0 #use to slice Ih table.
    h_index_matrix = sparse.matrix(n_refl, group_boundaries[1] - group_boundaries[0])

    for refl_idx, asu_idx in enumerate(asu_miller_index):
      if asu_idx != previous:
        total_refl_group_idx += 1
        refl_group_idx += 1
      if total_refl_group_idx == group_boundaries[block_idx + 1]:
        #save the previous group
        sub_refl_table = reflection_table[prev_refl_idx:refl_idx]
        sub_nz_weights = nonzero_weights[prev_refl_idx:refl_idx]
        sel1 = self._permutation_selection >= prev_refl_idx
        sel2 = self._permutation_selection < refl_idx
        block_sel = sel1 & sel2
        permuted_isel = self._permutation_selection.select(block_sel)
        selector = flex.bool(block_sel.size(), False)
        selector.set_selected(permuted_isel, True)
        col_selector = flex.size_t(range(0, h_index_matrix.non_zeroes))
        h_idx_T = h_index_matrix.transpose()
        reduced_h_idx_T = h_idx_T.select_columns(col_selector)
        reduced_h_idx = reduced_h_idx_T.transpose()
        self.add_Ihtable_block(selector, sub_refl_table, reduced_h_idx,
          sub_nz_weights)
        #start the next h_idx_matrix
        block_idx += 1
        h_index_matrix = sparse.matrix(n_refl, (group_boundaries[block_idx+1]
          - group_boundaries[block_idx]))
        prev_refl_idx = refl_idx
        refl_in_group_idx = 0
        refl_group_idx = 0
      h_index_matrix[refl_in_group_idx, refl_group_idx] = 1.0
      previous = asu_idx
      refl_in_group_idx += 1
    sub_refl_table = reflection_table[prev_refl_idx:]
    sub_nz_weights = nonzero_weights[prev_refl_idx:]
    block_sel = self._permutation_selection >= prev_refl_idx
    permuted_isel = self._permutation_selection.select(block_sel)
    selector = flex.bool(block_sel.size(), False)
    selector.set_selected(permuted_isel, True)
    col_selector = flex.size_t(range(0, h_index_matrix.non_zeroes))
    h_idx_T = h_index_matrix.transpose()
    reduced_h_idx_T = h_idx_T.select_columns(col_selector)
    reduced_h_idx = reduced_h_idx_T.transpose()
    self.add_Ihtable_block(selector, sub_refl_table, reduced_h_idx,
      sub_nz_weights)

    '''if self.free_set_sel:
      isel = (~self.free_set_sel).iselection()
      total_h_idx_matrix_T = total_h_idx_matrix.transpose()
      total_h_idx_matrix_T = total_h_idx_matrix_T.select_columns(isel)
      total_h_idx_matrix = total_h_idx_matrix_T.transpose()
      good_cols = flex.size_t([])
      for i, col in enumerate(total_h_idx_matrix.cols()):
        if col.non_zeroes != 0:
          good_cols.append(i)
      total_h_idx_matrix = total_h_idx_matrix.select_columns(good_cols)
    total_h_expand_matrix = total_h_idx_matrix.transpose()
    return total_h_idx_matrix, total_h_expand_matrix'''

  def add_Ihtable_block(self, selection, Ih_table, h_index_matrix,
    nonzero_weights_isel):
    """Add a new Ih table block to the list."""
    self._blocked_data_list.append(SingleIhTable(Ih_table, h_index_matrix,
      self.space_group, self.weighting_scheme, nonzero_weights_isel))
    self._blocked_selection_list.append(selection)

  def apply_selection_to_selection_list(self, i, sel):
    """"Apply a selection to the block selection list. Called by TargetScaler
    to reduce the data after matching equivalent reflections."""
    start_idx = 0
    for j, selection_list in self.blocked_selection_list.iteritems():
      current_sel_list = selection_list[i]
      n = len(current_sel_list)
      end_idx = start_idx + n
      new_sel_list = current_sel_list.select(sel[start_idx:end_idx])
      self.blocked_selection_list[j][i] = new_sel_list
      start_idx += n

  def sort_by_dataset_id(self):
    """Iterate through the blocks, sorting them to be ordered first by dataset
    id, then by asu_miller_index, rather than by just by asu_miller_index."""
    for block in self.blocked_data_list:
      if 'dataset_id' in block.Ih_table:
        #First sort the reflection table and reorder the h_index matrices.
        dataset_sort = flex.sort_permutation(block.Ih_table['dataset_id'])
        block.dataset_sort = dataset_sort
        block.Ih_table = block.Ih_table.select(dataset_sort)
        h_idx_T = block.h_index_matrix.transpose()
        per_h_idx_T = h_idx_T.select_columns(dataset_sort)
        block.h_index_matrix = per_h_idx_T.transpose()
        block.h_expand_matrix = block.h_index_matrix.transpose()
    #Now rearrange the selection lists.
    block_selections = {i:[] for i in range(self._n_datasets)}
    for sel_list, block in zip(self.blocked_selection_list,
      self.blocked_data_list):
      a = self._permutation_selection.select(sel_list)
      b = a.select(block.dataset_sort)
      for datasetid_ in range(self._n_datasets):
        sel = block.Ih_table['dataset_id'] == datasetid_
        block_selections[datasetid_].append(b.select(sel))
    cumulative_size = 0
    # Shift the indices so that the selection lists relate to the individual
    # datasets, rather than the order in the total dataset.
    for i in range(self._n_datasets):
      for j, block_sel in enumerate(block_selections[i]):
        block_selections[i][j] = block_sel - flex.size_t(block_sel.size(),
          cumulative_size)
      cumulative_size += sum(k.size() for k in block_selections[i])
    self._blocked_selection_list = block_selections

  def set_derivatives(self, derivatives):
    """Set the derivatives for each block."""
    for deriv, block in zip(derivatives, self.blocked_data_list):
      block.derivatives = deriv

  def set_inverse_scale_factors(self, new_scales):
    """Set the inverse scale factors for each block."""
    for scales, block in zip(new_scales, self.blocked_data_list):
      block.inverse_scale_factors = scales

  def calc_Ih(self):
    """Calculate the latest value of Ih in each block."""
    for block in self._blocked_data_list:
      block.calc_Ih()

  def update_weights(self):
    """Update weights in each block."""
    for block in self.blocked_data_list:
      block.update_weights()

  def update_weights_from_error_models(self):
    """Reset weights, to be used after the weights of individual Ih tables
    have been modified."""
    for block in self.blocked_data_list:
      block.update_weights_from_error_models()


'''def split_into_free_work(self, percentage, offset=0):
  """Split the dataset into a working and free set. The reflections chosen are
  a percentage of the symmetry equivalent groups, chosen in even steps starting
  at the offset, default 0."""
  n_groups = self.h_index_matrix.n_cols
  group_selection = [i for i in range(offset, n_groups, int(100/percentage))]
  free_set = flex.double(self.size, 0.0)
  for g in group_selection:
    free_set += self.h_index_matrix.col(g).as_dense_vector()
  free_set_sel = free_set != 0.0
  #select out a reflection table, and use this to initialise a new Ih_table.
  free_Ih_table = self._Ih_table.select(free_set_sel)
  free_Ih_table.set_flags(flex.bool(free_Ih_table.size(), False),
    free_Ih_table.flags.bad_for_scaling)
  #self._free_Ih_table = SingleIhTable(free_Ih_table, self.space_group)
  self._free_set_sel = free_set_sel
  self = self.select(~free_set_sel)
  msg = (('The dataset has been split into a working and free set. \n'
    'The free set contains {0} reflections, in {1} unique groups,\n'
    'and the working set contains {2} reflections in {3} unique groups.\n'
    ).format(self.free_Ih_table.size, self.free_Ih_table.h_index_matrix.n_cols,
  self.size, self.h_index_matrix.n_cols))
  logger.info(msg)'''
