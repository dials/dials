"""
Module of outlier rejection algorithms.
"""
import logging
from libtbx.utils import Sorry
from dials.array_family import flex
from dials.algorithms.scaling.Ih_table import IhTable

logger = logging.getLogger('dials')

def reject_outliers(reflection_tables, space_group, method='standard', zmax=9.0,
  target=None):
  """Helper function to act as interface to outlier algorithms."""
  if method == 'standard':
    refls = NormDevOutlierRejection(reflection_tables, space_group,
      zmax).return_reflection_tables()
  elif method == 'simple':
    refls = SimpleNormDevOutlierRejection(reflection_tables, space_group,
      zmax).return_reflection_tables()
  elif method == 'target':
    assert target is not None
    refls = TargetedOutlierRejection(reflection_tables, space_group,
      zmax, target).return_reflection_tables()
  else:
    raise Sorry("Invalid choice of outlier rejection method.")
  if len(reflection_tables) > 1:
    msg = ('Combined outlier rejection has been performed across multiple datasets, \n')
  else:
    msg = ('A round of outlier rejection has been performed, \n')
  n_outliers = sum([refl.get_flags(refl.flags.outlier_in_scaling).count(True)
    for refl in reflection_tables])
  msg += ('{0} outliers have been identified. \n'.format(n_outliers))
  logger.info(msg)
  return refls

class OutlierRejectionBase(object):
  """Base class for outlier rejection algorithms based on the use of the
  Ih_table datastructure."""

  def __init__(self, reflection_tables, space_group, zmax):
    self.reflection_tables = reflection_tables
    self.space_group = space_group
    self.zmax = zmax
    self.outliers_lists = [flex.size_t([]) for _ in self.reflection_tables]
    # First all outliers are unset, as the Ih_table will not
    # usee reflections flagged as outliers.
    self.unset_outlier_flags()
    self.do_outlier_rejection()
    self.set_outlier_flags()

  def return_reflection_tables(self):
    """Return the reflection table."""
    return self.reflection_tables

  def unset_outlier_flags(self):
    """Unset the outlier_in_scaling flag in the reflection table."""
    for refl_table in self.reflection_tables:
      already_outliers = refl_table.get_flags(
        refl_table.flags.outlier_in_scaling)
      already_sel = already_outliers == True
      if already_sel.count(True):
        refl_table.unset_flags(already_sel,
          refl_table.flags.outlier_in_scaling)

  def set_outlier_flags(self):
    """Set the outlier_in_scaling flag in the reflection table."""
    for i, refl_table in enumerate(self.reflection_tables):
      outlier_mask = flex.bool(refl_table.size(), False)
      outlier_mask.set_selected(self.outliers_lists[i], True)
      refl_table.set_flags(outlier_mask,
        refl_table.flags.outlier_in_scaling)

class TargetedOutlierRejection(OutlierRejectionBase):
  """Outlier rejection routine with a target of 'ideal' intensities."""

  def __init__(self, reflection_tables, space_group, zmax, targets):
    self.target_tables = targets
    super(TargetedOutlierRejection, self).__init__(
      reflection_tables, space_group, zmax)

  def do_outlier_rejection(self):
    """Do the outlier rejection algorithm."""
    outlier_indices_lists = self.round_of_outlier_rejection(
      self.reflection_tables, self.target_tables)
    for i, outlier_list in enumerate(outlier_indices_lists):
      self.outliers_lists[i].extend(outlier_list)

  def round_of_outlier_rejection(self, reflection_tables, targets):
    """One round of outlier rejection for all data in reflection table."""
    Ih_table = IhTable([(i, None) for i in reflection_tables], self.space_group,
      n_blocks=1).blocked_data_list[0]
    target_Ih_table = IhTable([(i, None) for i in targets], self.space_group,
      n_blocks=1).blocked_data_list[0]
    target_asu_Ih_dict = dict(zip(target_Ih_table.asu_miller_index,
      zip(target_Ih_table.Ih_values, target_Ih_table.variances)))
    Ih_table.Ih_table['target_Ih_value'] = flex.double(Ih_table.size, 0.0)
    Ih_table.Ih_table['target_Ih_sigmasq'] = flex.double(Ih_table.size, 0.0)
    for j, miller_idx in enumerate(Ih_table.asu_miller_index):
      if miller_idx in target_asu_Ih_dict:
        Ih_table.Ih_table['target_Ih_value'][j] = target_asu_Ih_dict[miller_idx][0]
        Ih_table.Ih_table['target_Ih_sigmasq'][j] = target_asu_Ih_dict[miller_idx][1]
    Ih_table.select(Ih_table.Ih_table['target_Ih_value'] != 0.0)
    norm_dev = (Ih_table.intensities - (
      Ih_table.inverse_scale_factors * Ih_table.Ih_table['target_Ih_value']))/ \
      ((Ih_table.variances + ((Ih_table.inverse_scale_factors**2) * \
      Ih_table.Ih_table['target_Ih_sigmasq']))**0.5)
    outliers_sel = flex.abs(norm_dev) > self.zmax

    nz_weights_isel = Ih_table.nonzero_weights#.iselection()
    outlier_indices = nz_weights_isel.select(outliers_sel)
    ids = Ih_table.Ih_table['dataset_id'].select(outliers_sel)

    outliers_list = []
    for i in range(len(reflection_tables)):
      sel = ids == i
      outlier_indices_for_i = outlier_indices.select(sel)
      outliers_list.append(outlier_indices_for_i)
    return outliers_list


class SimpleNormDevOutlierRejection(OutlierRejectionBase):
  """Outlier rejection algorithm using normalised deviations."""

  def do_outlier_rejection(self):
    """Do the outlier rejection algorithm."""
    outlier_indices_lists = self.round_of_outlier_rejection(
      self.reflection_tables)
    for i, outlier_list in enumerate(outlier_indices_lists):
      self.outliers_lists[i].extend(outlier_list)

  def round_of_outlier_rejection(self, reflection_tables):
    """One round of outlier rejection for all data in reflection table."""
    Ih_table = IhTable([(i, None) for i in reflection_tables], self.space_group,
      n_blocks=1).blocked_data_list[0]
    I = Ih_table.intensities
    g = Ih_table.inverse_scale_factors
    w = Ih_table.weights
    wgIsum = ((w * g * I) * Ih_table.h_index_matrix) * Ih_table.h_expand_matrix
    wg2sum = ((w * g * g) * Ih_table.h_index_matrix) * Ih_table.h_expand_matrix
    norm_dev = (I - (g * wgIsum/wg2sum))/(((1.0/w)+((g/wg2sum)**2))**0.5)
    outliers_sel = flex.abs(norm_dev) > self.zmax

    nz_weights_isel = Ih_table.nonzero_weights#.iselection()
    outlier_indices = nz_weights_isel.select(outliers_sel)
    ids = Ih_table.Ih_table['dataset_id'].select(outliers_sel)

    outliers_list = []
    for i in range(len(reflection_tables)):
      sel = ids == i
      outlier_indices_for_i = outlier_indices.select(sel)
      outliers_list.append(outlier_indices_for_i)
    return outliers_list


class NormDevOutlierRejection(OutlierRejectionBase):
  """Outlier rejection algorithm using normalised deviations."""

  def do_outlier_rejection(self):
    """Do the outlier rejection algorithm."""
    #First round of rejection
    outlier_indices_lists, other_potential_outliers_list = \
      self.round_of_outlier_rejection(self.reflection_tables)
    for i, outlier_list in enumerate(outlier_indices_lists):
      self.outliers_lists[i].extend(outlier_list)
    self.check_for_more_outliers(self.reflection_tables,
      other_potential_outliers_list)

  def check_for_more_outliers(self, reflection_tables,
    other_potential_outliers_list):
    """Recursive check for further outliers."""
    if not all([i.size() == 0 for i in other_potential_outliers_list]):
      # Select only the data where we want to look for further outliers.
      new_rts = []
      for i, refl in enumerate(reflection_tables):
        new_rts.append(refl.select(other_potential_outliers_list[i]))
      internal_outlier_indices, internal_other_potential_outliers = (
        self.round_of_outlier_rejection(new_rts))
      # Now the new outlier indices are relative to the new input,
      # need to transform back to original indices.
      new_other_potential_outliers_list = []
      for i, outl_list in enumerate(other_potential_outliers_list):
        new_outlier_indices = outl_list.select(
          internal_outlier_indices[i])
        new_other_potential_outliers_list.append(outl_list.select(
          internal_other_potential_outliers[0]))
        self.outliers_lists[i].extend(new_outlier_indices)
      self.check_for_more_outliers(reflection_tables,
        new_other_potential_outliers_list)

  def round_of_outlier_rejection(self, reflection_tables):
    """One round of outlier rejection for all data in reflection table."""
    Ih_table = IhTable([(i, None) for i in reflection_tables],
      self.space_group, n_blocks=1).blocked_data_list[0]
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
    sel = nh > 2 #could be > 1 if we want to calculate z_score for groups of 2
    wg2sum_others_sel = wg2sum_others.select(sel)
    wgIsum_others_sel = wgIsum_others.select(sel)
    g_sel = g.select(sel)
    I_sel = I.select(sel)
    w_sel = w.select(sel)

    norm_dev = (I_sel - (g_sel * wgIsum_others_sel/wg2sum_others_sel))/(
      ((1.0/w_sel)+((g_sel/wg2sum_others_sel)**2))**0.5)
    z_score = flex.abs(norm_dev)
    # Want an array same size as Ih table.
    all_z_scores = flex.double(Ih_table.size, 0.0)
    all_z_scores.set_selected(sel.iselection(), z_score)

    outliers = flex.bool(sel.size(), False)
    other_potential_outliers = flex.size_t([])
    for col in Ih_table.h_index_matrix.cols():
      if col.non_zeroes > 2:
        sel = col.as_dense_vector() > 0
        indices_of_group = sel.iselection()
        z_scores = all_z_scores.select(indices_of_group)
        max_z = max(z_scores)
        if max_z > self.zmax:
          max_selecter = z_scores == max_z
          outlier_index = indices_of_group.select(max_selecter)
          outliers.set_selected(outlier_index, True)
          other_indices = indices_of_group.select(~max_selecter)
          other_potential_outliers.extend(other_indices)

    #Now determine location of outliers w.r.t initial reflection table order.
    nz_weights_isel = Ih_table.nonzero_weights#.iselection()
    outlier_indices = nz_weights_isel.select(outliers)
    other_potential_outliers_indices = nz_weights_isel.select(
      other_potential_outliers)

    outlier_ids = Ih_table.Ih_table['dataset_id'].select(outliers)
    other_potential_ids = Ih_table.Ih_table['dataset_id'].select(
      other_potential_outliers)

    outliers_list = []
    other_potential_list = []
    for i in range(len(reflection_tables)):
      #sel = outlier_ids == i
      outliers_list.append(outlier_indices.select(outlier_ids == i))
      other_potential_list.append(other_potential_outliers_indices.select(
        other_potential_ids == i))
    return outliers_list, other_potential_list
