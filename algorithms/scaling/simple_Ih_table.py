from libtbx.containers import OrderedSet
from dials.array_family import flex
from cctbx import miller, crystal
from scitbx import sparse

def map_indices_to_asu(miller_indices, space_group):
  """Map the indices to the asymmetric unit."""
  crystal_symmetry = crystal.symmetry(space_group=space_group)
  miller_set = miller.set(crystal_symmetry=crystal_symmetry,
    indices=miller_indices, anomalous_flag=False)
  miller_set_in_asu = miller_set.map_to_asu()
  return miller_set_in_asu.indices()

def get_sorted_asu_indices(asu_indices, space_group):
  """Return the sorted asu indices and the permutation selection."""
  crystal_symmetry = crystal.symmetry(space_group=space_group)
  miller_set_in_asu = miller.set(crystal_symmetry=crystal_symmetry,
    indices=asu_indices, anomalous_flag=False)
  permuted = miller_set_in_asu.sort_permutation(by_value='packed_indices')
  sorted_asu_miller_index = asu_indices.select(permuted)
  return sorted_asu_miller_index, permuted

class simple_Ih_table(object):

  """
  The main datastructure used by scaling algorithms for performing operations
  on groups of symmetry equivalent reflections.

  The idea here is to split the data into blocks to allow parallelized
  computations, but within the blocks the data are sorted by dataset.
  In each block, there exists a block_selection_list which contains the indices
  for each dataset from the input reflection table.

  This class acts as a 'master' to setup the block structure and control access
  to the underlying blocks - only metadata is kept in this class after
  initialisation, the reflections etc are all contained in the blocks.
  To set data in the blocks, methods are provided by the master, e.g
  set_intensities(intensities, block_id) which are delegated down to the
  appropriate block.

  Attributes:
      space_group: The space group for the dataset.
      Ih_table_blocks (list): A list of Ih_table_block instances. All symmetry
          equivalent reflections are recorded in the same block, to allow
          splitting of the dataset for parallelized computations.
      nblocks (int): The number of blocks in the Ih_table_blocks list.
      blocked_selection_list (list): A list of lists. bsl[i][j] is the selection
         list for block i, dataset j.
      n_datasets: The number of input reflection tables used to make the Ih_table.
      size: The number of reflections across all blocks
      asu_index_dict (dict): A dictionary, key: asu_miller_index, value tuple
          containing group_id and block_id (where group id is the group index
          within its block).

  """

  id_ = "IhTable"

  def __init__(self, list_of_reflections_and_selections, space_group, nblocks=1):
    self.asu_index_dict = {}
    self.space_group = space_group
    self.nblocks = nblocks
    self.n_datasets = len(list_of_reflections_and_selections)
    self.Ih_table_blocks = []
    self.size = None
    self.properties_dict = {
      'n_unique_in_each_block' : [],
      'n_reflections_in_each_block' : {},
      'miller_index_boundaries' : []}
    self._determine_required_block_structures(
      list_of_reflections_and_selections, nblocks)
    self._create_empty_Ih_table_blocks()
    for i, dataset in enumerate(list_of_reflections_and_selections):
      self._add_dataset_to_blocks(i, dataset)
    for block in self.Ih_table_blocks:
      block.finished()
    self.blocked_selection_list = [block.block_selections for block in self.Ih_table_blocks]
    self.free_Ih_table = None

  def update_weights(self, block_id=None):
    pass

  def update_error_model(self, error_model):
    """Update the error model in the blocks."""
    for block in self.Ih_table_blocks:
      block.update_error_model(error_model)

  @property
  def blocked_data_list(self):
    return self.Ih_table_blocks

  def set_intensities(self, intensities, block_id):
    self.Ih_table_blocks[block_id].Ih_table['intensity'] = intensities

  def set_derivatives(self, derivatives, block_id):
    """Set the derivatives for each block."""
    self.Ih_table_blocks[block_id].derivatives = derivatives

  def set_inverse_scale_factors(self, new_scales, block_id):
    """Set the inverse scale factors for each block."""
    self.Ih_table_blocks[block_id].inverse_scale_factors = new_scales

  def set_variances(self, new_variances, block_id):
    """Set the inverse scale factors for each block."""
    self.Ih_table_blocks[block_id].Ih_table['variance'] = new_variances
    self.Ih_table_blocks[block_id].Ih_table['weights'] = 1.0 / new_variances

  def calc_Ih(self, block_id=None):
    """Calculate the latest value of Ih in each block."""
    if block_id:
      self.Ih_table_blocks[block_id].calc_Ih()
    else:
      for block in self.Ih_table_blocks:
        block.calc_Ih()

  def _determine_required_block_structures(self, list_of_reflections_and_selections, nblocks=1):
    """Extract the asu miller indices from the reflection table and
    add data to the asu_index_dict and properties dict."""
    joint_asu_indices = flex.miller_index()
    for refl_and_sel in list_of_reflections_and_selections:
      if not 'asu_miller_index' in refl_and_sel[0]:
        refl_and_sel[0]['asu_miller_index'] = map_indices_to_asu(
          refl_and_sel[0]['miller_index'], self.space_group)
      if refl_and_sel[1]:
        joint_asu_indices.extend(refl_and_sel[0]['asu_miller_index'].select(
          refl_and_sel[1]))
      else:
        joint_asu_indices.extend(refl_and_sel[0]['asu_miller_index'])
    sorted_joint_asu_indices, _ = get_sorted_asu_indices(
      joint_asu_indices, self.space_group)
    self.size = sorted_joint_asu_indices.size()
    asu_index_set = OrderedSet(sorted_joint_asu_indices)
    n_unique_groups = len(asu_index_set)

    #also record how many unique groups go into each block
    group_boundaries = [int(i*n_unique_groups/nblocks) for i in range(nblocks)]
    group_boundaries.append(n_unique_groups)

    next_boundary = group_boundaries[1]
    block_id = 0
    group_id_in_block_i = 0
    for i, index in enumerate(asu_index_set):
      if i == next_boundary:
        self.properties_dict['n_unique_in_each_block'].append(group_id_in_block_i)
        self.properties_dict['miller_index_boundaries'].append(index)
        block_id += 1
        next_boundary = group_boundaries[block_id+1]
        group_id_in_block_i = 0
      self.asu_index_dict[index] = (group_id_in_block_i, block_id)
      group_id_in_block_i += 1
    #record the number in the last block
    self.properties_dict['n_unique_in_each_block'].append(group_id_in_block_i)
    self.properties_dict['miller_index_boundaries'].append((0, 0, 0)) # to avoid bounds checking when in last group
    #need to know how many reflections will be in each block also

    block_id = 0
    idx_prev = 0
    boundary = self.properties_dict['miller_index_boundaries'][0]
    for i, index in enumerate(sorted_joint_asu_indices):
      if index == boundary:
        n_in_prev_group = i - idx_prev
        self.properties_dict['n_reflections_in_each_block'][block_id] = n_in_prev_group
        block_id += 1
        boundary = self.properties_dict['miller_index_boundaries'][block_id]
        idx_prev = i
    self.properties_dict['n_reflections_in_each_block'][block_id] = \
      len(sorted_joint_asu_indices) - idx_prev

  def _create_empty_Ih_table_blocks(self):
    for n in range(self.nblocks):
      n_refl_in_block = self.properties_dict['n_reflections_in_each_block'][n]
      n_groups_in_block = self.properties_dict['n_unique_in_each_block'][n]
      self.Ih_table_blocks.append(Ih_table_block(n_groups=n_groups_in_block,
        n_refl=n_refl_in_block, ndatasets=self.n_datasets))

  def _add_dataset_to_blocks(self, dataset_id, reflections_and_selection):
    selection = reflections_and_selection[1]
    reflections = reflections_and_selection[0]
    if selection:
      asu_indices = reflections['asu_miller_index'].select(selection)
    else:
      asu_indices = reflections['asu_miller_index']
      selection = flex.bool(reflections.size(), True)
    sorted_asu_indices, perm = get_sorted_asu_indices(
      asu_indices, self.space_group)
    r = flex.reflection_table()
    r['intensity'] = reflections['intensity']
    r['asu_miller_index'] = reflections['asu_miller_index']
    r['variance'] = reflections['variance']
    r['inverse_scale_factor'] = reflections['inverse_scale_factor']
    if selection:
      r = r.select(selection).select(perm)
      r['loc_indices'] = selection.iselection().select(perm)
    else:
      n = r.size()
      r = r.select(perm)
      r['loc_indices'] = flex.size_t(range(n)).select(perm)
    r['dataset_id'] = flex.int(r.size(), dataset_id)
    # if data are sorted by asu_index, then up until boundary, should be in same
    # block (still need to read group_id though)

    #sort data, get group ids and block_ids
    group_ids = flex.int([])
    boundary = self.properties_dict['miller_index_boundaries'][0]
    boundary_id = 0
    boundaries_for_this_datset = [0]#use to slice
    # make this a c++ method for speed?
    for i, index in enumerate(sorted_asu_indices):
      if index == boundary:
        boundaries_for_this_datset.append(i)
        boundary_id += 1
        boundary = self.properties_dict['miller_index_boundaries'][boundary_id]
      group_id, _ = self.asu_index_dict[index]
      group_ids.append(group_id)
    boundaries_for_this_datset.append(len(sorted_asu_indices))

    ##FIXME make add_data a c++ method as constructing h_index matrix
    # so now have group ids as well for individual dataset
    for i, val in enumerate(boundaries_for_this_datset[:-1]):
      start = val
      end = boundaries_for_this_datset[i+1]
      self.Ih_table_blocks[i].add_data(dataset_id, group_ids[start:end], r[start:end])
      # add data to block block_id, add a line to h_index_matix of that block in loc group_id


class Ih_table_block(object):

  def __init__(self, n_groups, n_refl, ndatasets=1):
    self.Ih_table = flex.reflection_table()
    self.block_selections = [None] * ndatasets
    self.h_index_matrix = sparse.matrix(n_refl, n_groups)
    self.next_row = 0
    self.next_dataset = 0
    self._n_h = None
    self.h_expand_matrix = None
    self.derivatives = None

  def add_data(self, dataset_id, group_ids, reflections):
    #make h_index matrix
    assert dataset_id == self.next_dataset
    for i, id_ in enumerate(group_ids):
      self.h_index_matrix[i+self.next_row, id_] = 1.0
    self.next_row += len(group_ids)
    self.next_dataset += 1
    self.Ih_table.extend(reflections)
    self.block_selections[dataset_id] = reflections['loc_indices']

  def select(self, sel):
    """Select a subset of the data, returning a new Ih_table_block object."""
    Ih_table = self.Ih_table.select(sel)
    h_idx_sel = self.h_expand_matrix.select_columns(sel.iselection())
    reduced_h_idx = h_idx_sel.transpose()
    unity = flex.double(reduced_h_idx.n_rows, 1.0)
    nz_col_sel = (unity * reduced_h_idx) > 0
    h_index_matrix = reduced_h_idx.select_columns(nz_col_sel.iselection())
    h_expand = h_index_matrix.transpose()
    newtable = Ih_table_block(n_groups=0, n_refl=0)
    newtable.Ih_table = Ih_table
    newtable.h_expand_matrix = h_expand
    newtable.h_index_matrix = h_index_matrix
    return newtable

  def finished(self):
    self.h_index_matrix.compact()
    self.h_expand_matrix = self.h_index_matrix.transpose()
    self.Ih_table['weights'] = 1.0/self.Ih_table['variance']

  def calc_Ih(self):
    """Calculate the current best estimate for I for each reflection group."""
    scale_factors = self.Ih_table['inverse_scale_factor']
    gsq = ((scale_factors**2) * self.Ih_table['weights'])
    sumgsq = gsq * self.h_index_matrix
    gI = ((scale_factors * self.Ih_table['intensity']) * self.Ih_table['weights'])
    sumgI = gI * self.h_index_matrix
    Ih = sumgI/sumgsq
    self.Ih_table['Ih_values'] = Ih * self.h_expand_matrix

  @property
  def inverse_scale_factors(self):
    """"The inverse scale factors of the reflections."""
    return self.Ih_table['inverse_scale_factor']

  @inverse_scale_factors.setter
  def inverse_scale_factors(self, new_scales):
    if new_scales.size() != self.size:
      assert 0, """attempting to set a new set of scale factors of different
      length than previous assignment: was %s, attempting %s""" % (
        self.inverse_scale_factors.size(), new_scales.size())
    else:
      self.Ih_table['inverse_scale_factor'] = new_scales

  def update_error_model(self, error_model):
    """Update the scaling weights based on an error model."""
    sigmaprimesq = error_model.update_variances(self.Ih_table['variance'], self.Ih_table['intensity'])
    self.Ih_table['weights'] = 1.0/sigmaprimesq

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

  @property
  def variances(self):
    """The initial variances of the reflections."""
    return self.Ih_table['variance']

  @property
  def intensities(self):
    """The unscaled reflection intensities."""
    return self.Ih_table['intensity']

  @property
  def Ih_values(self):
    """The estimated intensities of symmetry equivalent reflections."""
    return self.Ih_table['Ih_values']

  @property
  def weights(self):
    """The weights that will be used in scaling."""
    return self.Ih_table['weights']

  @weights.setter
  def weights(self, new_weights):
    if new_weights.size() != self.size:
      assert 0, '''attempting to set a new set of weights of different
      length than previous assignment: was %s, attempting %s''' % (
        self.size, new_weights.size())
    self.Ih_table['weights'] = new_weights

  @property
  def size(self):
    """Return the length of the stored Ih_table (a reflection table)."""
    return self.Ih_table.size()

  @property
  def asu_miller_index(self):
    """Return the length of the stored Ih_table (a reflection table)."""
    return self.Ih_table['asu_miller_index']
