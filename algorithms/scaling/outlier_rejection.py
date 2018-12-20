"""
Module of outlier rejection algorithms.
"""
import abc
import logging
from libtbx.utils import Sorry
from scitbx.array_family import flex
from dials.algorithms.scaling.simple_Ih_table import simple_Ih_table as IhTable
from dials_scaling_ext import determine_outlier_indices

logger = logging.getLogger('dials')

def reject_outliers(reflection_table, experiment, method='standard', zmax=6.0):
  """Run an outlier algorithm on symmetry-equivalent intensities.

  This method runs an intensity-based outlier rejection algorithm, comparing
  the deviations from the weighted mean in groups of symmetry equivalent
  reflections. The outliers are determined and the outlier_in_scaling flag
  is set in the reflection table.

  The values 'intensity' and 'variance' must be set in the reflection table;
  these should be corrected but unscaled values, as an 'inverse_scale_factor'
  will be applied during outlier rejection if this is present in the reflection
  table. The reflection table should also be prefiltered (e.g. not-integrated
  reflections should not be present) as no further filtering is done on the
  input table.

  Args:
      reflection_table: A reflection table.
      experiment: A single experiment object.
      method (str): Name (alias) of outlier rejection algorithm to use.
      zmax (float): Normalised deviation threshold for classifying an outlier.

  Returns:
      reflection_table: The input table with the outlier_in_scaling flag set.

  Raises:
      Sorry: if the reflection table does not contain 'intensity' and 'variance'.
  """
  if not 'intensity' in reflection_table or not 'variance' in reflection_table:
    raise Sorry("""The reflection table does not contain columns with the names
      'intensity' and 'variance'""")

  if not 'inverse_scale_factor' in reflection_table:
    reflection_table['inverse_scale_factor'] = flex.double(
      reflection_table.size(), 1.0)

  Ih_table = IhTable([reflection_table],
    experiment.crystal.get_space_group(), nblocks=1)
  outlier_indices = determine_outlier_index_arrays(
    Ih_table, method=method, zmax=zmax)[0]

  # Unset any existing outlier flags before setting the new ones
  reflection_table.unset_flags(reflection_table.get_flags(
    reflection_table.flags.outlier_in_scaling),
    reflection_table.flags.outlier_in_scaling)
  reflection_table.set_flags(outlier_indices,
    reflection_table.flags.outlier_in_scaling)

  return reflection_table


def determine_outlier_index_arrays(Ih_table, method='standard', zmax=6.0,
  target=None):
  """Setup and run an outlier algorithm and return the outlier indices.

  Args:
      Ih_table: A dials.algorithms.scaling.simple_Ih_table.simple_Ih_table.
      method (str): Name (alias) of outlier rejection algorithm to use. If
          method='target', then the optional argument 'target' must also
          be specified.
      zmax (float): Normalised deviation threshold for classifying an outlier.
      target (Optional[Ih_table]): An Ih_table to use to obtain target Ih for
          outlier rejectiob, if method='target'.

  Returns:
      outlier_index_arrays (list): A list of flex.size_t arrays, with one
          array per dataset that was used to create the Ih_table. Importantly,
          the indices are the indices of the reflections in the initial
          reflection table used to create the Ih_table, not the indices of the
          data in the Ih_table.

  """
  if method == 'standard':
    outlier_index_arrays = NormDevOutlierRejection(Ih_table,
      zmax).final_outlier_arrays
  elif method == 'simple':
    outlier_index_arrays = SimpleNormDevOutlierRejection(Ih_table,
      zmax).final_outlier_arrays
  elif method == 'target':
    assert target is not None
    outlier_index_arrays = TargetedOutlierRejection(Ih_table,
      zmax, target).final_outlier_arrays
  elif method is None:
    return [flex.size_t([]) for _ in range(Ih_table.n_datasets)]
  else:
    raise Sorry("Invalid choice of outlier rejection method.")
  if Ih_table.n_datasets > 1:
    msg = ('Combined outlier rejection has been performed across multiple datasets, \n')
  else:
    msg = ('A round of outlier rejection has been performed, \n')
  n_outliers = sum([len(i) for i in outlier_index_arrays])
  msg += ('{0} outliers have been identified. \n'.format(n_outliers))
  logger.info(msg)
  return outlier_index_arrays


class OutlierRejectionBase(object):
  """
  Base class for outlier rejection algorithms based on the use of the
  Ih_table datastructure.

  Subclasses must define the do_outlier_rejection method, which must
  add the indices of outliers to the outlier_indices attribute.

  Attributes:
      Ih_table_block: An Ih_table_block of reflections to test for outliers.
          This may be reduced in size during the algorithm.
      n_datasets (int): The number of reflection tables used to create the
          Ih_table.
      block_selections (list): A list of flex.size_t arrays that relate the
          order of the Ih_table data to the initial reflection tables.
      ids: A flex.int array of dataset ids (0..n-1).
      zmax (float): Normalised deviation threshold for classifying an outlier.
      outlier_indices: A flex.size_t array of outlier indices w.r.t. the
          Ih_table data order.
      final_outlier_arrays (list): A list of flex.size_t arrays of outlier
          indices w.r.t. the order of the initial reflection tables used to
          create the Ih_table.

  """

  __metaclass__ = abc.ABCMeta

  def __init__(self, Ih_table, zmax):
    assert Ih_table.n_work_blocks == 1, """
Outlier rejection algorithms require an Ih_table with nblocks = 1"""
    # Note: could be possible to code for nblocks > 1
    self.Ih_table_block = Ih_table.blocked_data_list[0]
    self.n_datasets = Ih_table.n_datasets
    self.block_selections = Ih_table.blocked_selection_list[0]
    self.ids = self.Ih_table_block.Ih_table['dataset_id']
    self.zmax = zmax
    self.outlier_indices = flex.size_t([])
    self.do_outlier_rejection()
    self.final_outlier_arrays = self.determine_outlier_indices()

  def determine_outlier_indices(self):
    """Transform the outlier indices w.r.t the Ih_table to outlier indices
    w.r.t the initial reflection tables used to create the Ih_table, separated
    by reflection table.

    Returns:
        final_outlier_arrays (list): A list of flex.size_t arrays of
            outlier indices w.r.t. the order of the data in the initial
            reflection tables used to create the Ih_table.

    """
    if self.n_datasets == 1:
      return [self.block_selections[0].select(self.outlier_indices)]
    final_outlier_arrays = []
    ids = self.ids.select(self.outlier_indices)
    offset = 0
    for i in range(self.n_datasets):
      outlier_array_i = self.outlier_indices.select(ids == i) - offset
      final_outlier_arrays.append(
        self.block_selections[i].select(outlier_array_i))
      offset += self.block_selections[i].size()
    return final_outlier_arrays

  @abc.abstractmethod
  def do_outlier_rejection(self):
    """Add indices (w.r.t. the Ih_table data) to self.outlier_indices"""


class TargetedOutlierRejection(OutlierRejectionBase):
  """Outlier rejection routine with a target of 'ideal' intensities."""

  def __init__(self, reflection_tables, zmax, target):
    assert target.n_work_blocks == 1, """
Targeted outlier rejection requires a target Ih_table with nblocks = 1"""
    self.target_Ih_table_block = target.blocked_data_list[0]
    self.target_Ih_table_block.calc_Ih()
    super(TargetedOutlierRejection, self).__init__(
      reflection_tables, zmax)

  def do_outlier_rejection(self):
    Ih_table = self.Ih_table_block
    target = self.target_Ih_table_block
    target_asu_Ih_dict = dict(zip(target.asu_miller_index,
      zip(target.Ih_values, target.variances)))
    Ih_table.Ih_table['target_Ih_value'] = flex.double(Ih_table.size, 0.0)
    Ih_table.Ih_table['target_Ih_sigmasq'] = flex.double(Ih_table.size, 0.0)
    for j, miller_idx in enumerate(Ih_table.asu_miller_index):
      if miller_idx in target_asu_Ih_dict:
        Ih_table.Ih_table['target_Ih_value'][j] = target_asu_Ih_dict[miller_idx][0]
        Ih_table.Ih_table['target_Ih_sigmasq'][j] = target_asu_Ih_dict[miller_idx][1]

    nz_sel = Ih_table.Ih_table['target_Ih_value'] != 0.0
    Ih_table = Ih_table.select(nz_sel)
    norm_dev = (Ih_table.intensities - (
      Ih_table.inverse_scale_factors * Ih_table.Ih_table['target_Ih_value']))/ \
      ((Ih_table.variances + ((Ih_table.inverse_scale_factors**2) * \
      Ih_table.Ih_table['target_Ih_sigmasq']))**0.5)
    outliers_sel = flex.abs(norm_dev) > self.zmax
    outliers_isel = nz_sel.iselection().select(outliers_sel)
    self.outlier_indices.extend(outliers_isel)


class SimpleNormDevOutlierRejection(OutlierRejectionBase):
  """
  Outlier rejection algorithm using normalised deviations from the
  weighted mean of all reflections in each group.
  """

  def do_outlier_rejection(self):
    Ih_table = self.Ih_table_block
    I = Ih_table.intensities
    g = Ih_table.inverse_scale_factors
    w = Ih_table.weights
    wgIsum = ((w * g * I) * Ih_table.h_index_matrix) * Ih_table.h_expand_matrix
    wg2sum = ((w * g * g) * Ih_table.h_index_matrix) * Ih_table.h_expand_matrix

    # guard against zero divison errors - can happen due to rounding errors
    # or bad data giving g values are very small
    zero_sel = (wg2sum == 0.0)
    # set as one for now, then mark as outlier below. This will only affect if
    # g is near zero, if w is zero then throw an assertionerror.
    wg2sum.set_selected(zero_sel, 1.0)

    assert w.all_gt(0) # guard against division by zero
    norm_dev = (I - (g * wgIsum/wg2sum))/(((1.0/w)+((g/wg2sum)**2))**0.5)
    norm_dev.set_selected(zero_sel, 1000) # to trigger rejection
    outliers_sel = flex.abs(norm_dev) > self.zmax

    self.outlier_indices.extend(outliers_sel.iselection())


class NormDevOutlierRejection(OutlierRejectionBase):
  """
  Outlier rejection algorithm using normalised deviations from the
  weighted mean of reflections in each group, excluding the test
  reflection.
  """

  def do_outlier_rejection(self):
    outlier_indices, other_potential_outliers = \
      self.round_of_outlier_rejection()
    self.outlier_indices.extend(outlier_indices)
    if other_potential_outliers:
      good_sel = flex.bool(self.Ih_table_block.Ih_table.size(), False)
      good_sel.set_selected(other_potential_outliers, True)
      self.Ih_table_block = self.Ih_table_block.select(good_sel)
      self.check_for_more_outliers(other_potential_outliers)

  def check_for_more_outliers(self, other_potential_outliers):
    """Recursive check for further outliers.

    Each iteration creates a new reduced-size Ih_table_block, which retains
    only symmetry groups that need further testing. Outlier indices must be
    transformed to give indices with respect to the initial Ih_table_block.

    Args:
        other_potential_outliers: A flex.size_t array of indices with respect
            to the initial Ih_table data

    """
    # Find outlier indices with respect to reduced Ih_table block
    internal_outlier_indices, internal_other_potential_outliers = (
      self.round_of_outlier_rejection())
    outliers_wrt_original = other_potential_outliers.select(
      internal_outlier_indices)
    self.outlier_indices.extend(outliers_wrt_original)
    new_other_potential_outliers = other_potential_outliers.select(
      internal_other_potential_outliers)# still wrt original Ih_table data

    if new_other_potential_outliers:
      good_sel = flex.bool(self.Ih_table_block.size, False)
      good_sel.set_selected(internal_other_potential_outliers, True)
      self.Ih_table_block = self.Ih_table_block.select(good_sel)
      self.check_for_more_outliers(new_other_potential_outliers)

  def round_of_outlier_rejection(self):
    """Calculate normal deviations from the data in the Ih_table

    Returns:
        (tuple): tuple containing:
            outlier_indices: A flex.size_t array of outlier indices w.r.t
                the current Ih_table
            other_potential_outliers: A flex.size_t array of indices from
                the symmetry groups where outliers were found, excluding the
                indices of the outliers themselves (indices w.r.t current
                Ih_table).

    """
    Ih_table = self.Ih_table_block
    I = Ih_table.intensities
    g = Ih_table.inverse_scale_factors
    w = Ih_table.weights
    wgIsum = ((w * g * I) * Ih_table.h_index_matrix) * Ih_table.h_expand_matrix
    wg2sum = ((w * g * g) * Ih_table.h_index_matrix) * Ih_table.h_expand_matrix
    wgIsum_others = wgIsum - (w * g * I)
    wg2sum_others = wg2sum - (w * g * g)
    # Now do the rejection analyis if n_in_group > 2
    Ih_table.calc_nh()
    nh = Ih_table.n_h
    sel = nh > 2
    wg2sum_others_sel = wg2sum_others.select(sel)
    wgIsum_others_sel = wgIsum_others.select(sel)

    # guard against zero divison errors - can happen due to rounding errors
    # or bad data giving g values are very small
    zero_sel = (wg2sum_others_sel == 0.0)
    # set as one for now, then mark as outlier below. This will only affect if
    # g is near zero, if w is zero then throw an assertionerror.
    wg2sum_others_sel.set_selected(zero_sel, 1.0)
    g_sel = g.select(sel)
    I_sel = I.select(sel)
    w_sel = w.select(sel)

    assert w_sel.all_gt(0) # guard against division by zero
    norm_dev = (I_sel - (g_sel * wgIsum_others_sel/wg2sum_others_sel))/(
      ((1.0/w_sel)+((g_sel/wg2sum_others_sel)**2))**0.5)
    norm_dev.set_selected(zero_sel, 1000) # to trigger rejection
    z_score = flex.abs(norm_dev)
    # Want an array same size as Ih table.
    all_z_scores = flex.double(Ih_table.size, 0.0)
    all_z_scores.set_selected(sel.iselection(), z_score)
    outlier_indices, other_potential_outliers = determine_outlier_indices(
      Ih_table.h_index_matrix, all_z_scores, self.zmax)
    return outlier_indices, other_potential_outliers
