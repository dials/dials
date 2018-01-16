from __future__ import print_function
from dials.array_family import flex
import numpy as np

import logging
logger = logging.getLogger('dials.scale')

class Weighting(object):
  ''' This class defines a weighting object that takes in a reflection table,
  gives initial weights of 1/variance and has methods to set the weights of
  certain reflections to zero.'''
  def __init__(self, reflection_table):
    '''set initial weighting to be a statistical weighting'''
    self._scale_weighting = 1.0/reflection_table['variance']

  @property
  def weights(self):
    return self._scale_weighting

  @weights.setter
  def weights(self, new_weights):
    self._scale_weighting = new_weights

  def tukey_biweighting(self, Ih_table):
    z_score = flex.double([])
    zmax = 6.0
    for i, _ in enumerate(Ih_table.h_index_counter_array):
      h_idx_cumul = Ih_table.h_index_cumulative_array[i:i+2]
      Ihls = Ih_table.intensities[h_idx_cumul[0]:h_idx_cumul[1]]
      #gs = Ih_table.inverse_scale_factors[h_idx_cumul[0]:h_idx_cumul[1]]
      #ws = Ih_table.weights[h_idx_cumul[0]:h_idx_cumul[1]]
      var = Ih_table.variances[h_idx_cumul[0]:h_idx_cumul[1]]
      med = np.median(Ihls)
      sigma = max([np.median(var**0.5), np.median(Ihls - med)])
      z = (Ihls - med) / sigma
      z_score.extend(z)
    self.weights = (1.0 - ((z_score/zmax)**2))**2
    sel = self.weights < 0.0
    self.weights.set_selected(sel, 0.0)

  def set_unity_weighting(self, reflection_table):
    '''method to weight each reflection equally'''
    self.weights = flex.double([1.0]*len(reflection_table['variance']))

  def apply_Isigma_cutoff(self, reflection_table, ratio):
    '''method to set a zero weight below an I/sigma cutoff'''
    Ioversigma = reflection_table['intensity']/(reflection_table['variance']**0.5)
    sel = Ioversigma <= ratio
    self.weights.set_selected(sel, 0.0)

  def apply_dmin_cutoff(self, reflection_table, d_cutoff_value):
    '''method to set a zero weight below an d-value cutoff'''
    sel = reflection_table['d'] <= d_cutoff_value
    self.weights.set_selected(sel, 0.0)

  def remove_wilson_outliers(self, reflection_table):
    '''method to set a zero weight for any outliers as determined by
    Wilson_outlier_test.py'''
    if 'wilson_outlier_flag' in reflection_table:
      sel = reflection_table['wilson_outlier_flag']
      self.weights.set_selected(sel, 0.0)

  def apply_aimless_error_model(self, reflection_table, error_params):
    '''applies scaling factors to the errors of the intensities'''
    msg = ('Applying an error model to the variances used for scaling {sep}'
      'with the error model parameters {0:.5f}, {1:.5f}. {sep}').format(error_params[0],
      error_params[1], sep='\n')
    logger.info(msg)
    sel = self.weights != 0.0
    nonzero_weights = self.weights.select(sel)
    nz_intensities = reflection_table.select(sel)['intensity']
    sigmaprime = error_params[0] * (((1.0/nonzero_weights)
                                     + ((error_params[1] * nz_intensities)**2))**0.5)
    new_weights = 1.0/(sigmaprime**2)
    self.weights.set_selected(sel, new_weights)
