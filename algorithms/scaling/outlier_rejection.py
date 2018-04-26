"""
Module of uoutlier rejection algorithms.
"""
from __future__ import print_function
from dials.array_family import flex
from dials.algorithms.scaling.Ih_table import SingleIhTable


def reject_outliers(reflection_table, space_group, params):
  """Helper function to act as interface to outlier algorithms."""
  if params.scaling_options.outlier_rejection == 'standard':
    refl = NormDevOutlierRejection(reflection_table, space_group,
      params.scaling_options.outlier_zmax).return_reflection_table()
    return refl
  elif params.scaling_options.outlier_rejection == 'simple':
    refl = SimpleNormDevOutlierRejection(reflection_table, space_group,
      params.scaling_options.outlier_zmax).return_reflection_table()
    return refl
  return reflection_table

class OutlierRejectionBase(object):
  """Base class for outlier rejection algorithms based on the use of the
  Ih_table datastructure."""

  def __init__(self, reflection_table, space_group, zmax):
    self.reflection_table = reflection_table
    self.space_group = space_group
    self.zmax = zmax
    self.outliers_list = flex.size_t([])
    # First all outliers are unset, as the Ih_table will not
    # usee reflections flagged as outliers.
    self.unset_outlier_flags()
    self.do_outlier_rejection()
    self.set_outlier_flags()

  def return_reflection_table(self):
    """Return the reflection table."""
    return self.reflection_table

  def unset_outlier_flags(self):
    """Unset the outlier_in_scaling flag in the reflection table."""
    already_outliers = self.reflection_table.get_flags(
      self.reflection_table.flags.outlier_in_scaling)
    already_sel = already_outliers == True
    if already_sel.count(True):
      self.reflection_table.unset_flags(already_sel,
        self.reflection_table.flags.outlier_in_scaling)

  def set_outlier_flags(self):
    """Set the outlier_in_scaling flag in the reflection table."""
    outlier_mask = flex.bool(self.reflection_table.size(), False)
    outlier_mask.set_selected(self.outliers_list, True)
    self.reflection_table.set_flags(outlier_mask,
      self.reflection_table.flags.outlier_in_scaling)


class SimpleNormDevOutlierRejection(OutlierRejectionBase):
  """Outlier rejection algorithm using normalised deviations."""

  def do_outlier_rejection(self):
    """Do the outlier rejection algorithm."""
    outlier_indices = self.round_of_outlier_rejection(self.reflection_table)
    self.outliers_list.extend(outlier_indices)

  def round_of_outlier_rejection(self, reflection_table):
    """One round of outlier rejection for all data in reflection table."""
    Ih_table = SingleIhTable(reflection_table, self.space_group)
    I = Ih_table.intensities
    g = Ih_table.inverse_scale_factors
    w = Ih_table.weights
    wgIsum = ((w * g * I) * Ih_table.h_index_matrix) * Ih_table.h_expand_matrix
    wg2sum = ((w * g * g) * Ih_table.h_index_matrix) * Ih_table.h_expand_matrix
    norm_dev = (I - (g * wgIsum/wg2sum))/(((1.0/w)+((g/wg2sum)**2))**0.5)
    z_score = (norm_dev**2)**0.5
    outliers_sel = z_score > self.zmax

    nz_weights_isel = Ih_table.nonzero_weights.iselection()
    outlier_indices = nz_weights_isel.select(outliers_sel)
    return outlier_indices


class NormDevOutlierRejection(OutlierRejectionBase):
  """Outlier rejection algorithm using normalised deviations."""

  def do_outlier_rejection(self):
    """Do the outlier rejection algorithm."""
    #First round of rejection
    outlier_indices, other_potential_outliers = self.round_of_outlier_rejection(
      self.reflection_table)
    self.outliers_list.extend(outlier_indices)
    self.check_for_more_outliers(self.reflection_table, other_potential_outliers)

  def check_for_more_outliers(self, reflection_table, other_potential_outliers):
    """Recursive check for further outliers."""
    if other_potential_outliers:
      # Select only the data where we want to look for further outliers.
      new_rt = reflection_table.select(other_potential_outliers)
      internal_outlier_indices, internal_other_potential_outliers = (
        self.round_of_outlier_rejection(new_rt))
      # Now the new outlier indices are relative to the new input,
      # need to transform back to original indices.
      new_outlier_indices = other_potential_outliers.select(
        internal_outlier_indices)
      new_other_potential_outliers = other_potential_outliers.select(
        internal_other_potential_outliers)
      self.outliers_list.extend(new_outlier_indices)
      self.check_for_more_outliers(reflection_table, new_other_potential_outliers)

  def round_of_outlier_rejection(self, reflection_table):
    """One round of outlier rejection for all data in reflection table."""
    Ih_table = SingleIhTable(reflection_table, self.space_group)
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
    z_score = (norm_dev**2)**0.5
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
    nz_weights_isel = Ih_table.nonzero_weights.iselection()
    outlier_indices = nz_weights_isel.select(outliers)
    other_potential_outliers_indices = nz_weights_isel.select(
      other_potential_outliers)

    return outlier_indices, other_potential_outliers_indices
