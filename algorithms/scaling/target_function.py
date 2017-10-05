'''
Classes that take in a data manager object and return a residual and gradient
vector for minimisation for different scaling parameterisations.
'''
import numpy as np
import copy
from dials.array_family import flex

class target_function(object):
  '''Superclass that takes in a data manager object and returns a residual
  and gradient function for minimisation. Gradient function must be
  overloaded by the subclass.'''
  def __init__(self, data_manager_object):
    self.data_manager = data_manager_object
    #self.parameters = parameters

  def calculate_residual(self):
    '''returns a residual vector'''
    intensities = self.data_manager.sorted_reflections['intensity']
    scale_factors = self.data_manager.sorted_reflections['inverse_scale_factor']
    Ih_values = self.data_manager.sorted_reflections['Ih_values']
    weights = self.data_manager.weights_for_scaling.get_weights()
    R = ((((intensities - (scale_factors * Ih_values))**2) * weights))
    return R

  def calculate_gradient(self):
    assert 0, 'parameterization-specific gradient function must be defined in a subclass'

  def return_targets(self):
    '''return residual and gradient arrays'''
    #self.data_manager.update_for_minimisation(self.parameters)
    return self.calculate_residual(), self.calculate_gradient()

class xds_target_function(target_function):
  '''Subclass that takes a data manager object and returns a residual and
  gradient function for a xds-like scaling parameterisation.'''

  def calculate_gradient(self):
    gradient = flex.double([])
    intensities = self.data_manager.sorted_reflections['intensity']
    scale_factors = self.data_manager.sorted_reflections['inverse_scale_factor']
    Ih_values = self.data_manager.sorted_reflections['Ih_values']
    scaleweights = self.data_manager.weights_for_scaling.get_weights()
    gsq = ((scale_factors)**2) * scaleweights#/ variances
    sumgsq = np.add.reduceat(gsq, self.data_manager.h_index_cumulative_array[:-1])
    rhl = intensities - (Ih_values * scale_factors)
    num = len(intensities)
    dIh = ((intensities * scaleweights) - (Ih_values * 2.0 * scale_factors * scaleweights))

    for i in range(len(self.data_manager.g_parameterisation[self.data_manager.active_parameterisation]['parameterisation'])):
      dIh_g = (dIh * self.data_manager.g_parameterisation[
               self.data_manager.active_parameterisation]['derivatives'][i*num:(i+1)*num])
      dIh_g = np.add.reduceat(dIh_g, self.data_manager.h_index_cumulative_array[:-1])/sumgsq
      dIh_g = flex.double(np.repeat(dIh_g, self.data_manager.h_index_counter_array))
      drdp = -((Ih_values * self.data_manager.g_parameterisation[
                  self.data_manager.active_parameterisation]['derivatives'][i*num:(i+1)*num])
               + (scale_factors * dIh_g))
      grad = (2.0 * rhl * scaleweights * drdp)
      gradient.append(flex.sum(grad))
    return gradient

  def calculate_gradient_v2(self):
    n_gbins = len(self.data_manager.g_decay)
    intensities = self.data_manager.sorted_reflections['intensity']
    num = len(intensities)
    scale_factors = self.data_manager.sorted_reflections['inverse_scale_factor']
    Ih_values = self.data_manager.sorted_reflections['Ih_values']
    scaleweights = self.data_manager.weights_for_scaling.get_weights()
    gsq = ((scale_factors)**2) * scaleweights#/ variances
    sumgsq = np.add.reduceat(gsq, self.data_manager.h_index_cumulative_array[:-1])
    sumgsq_tiled = flex.double(np.tile(sumgsq, n_gbins))
    rhl = intensities - (Ih_values * scale_factors)

    rhl_tiled = flex.double(np.tile(rhl, n_gbins))
    scaleweights_tiled = flex.double(np.tile(scaleweights, n_gbins))
    Ih_values_tiled = flex.double(np.tile(Ih_values, n_gbins))
    scale_factors_tiled = flex.double(np.tile(scale_factors, n_gbins))
    gradient_reducer = []
    for i in range(0,n_gbins):
      gradient_reducer.append(i*num)

    dIh = ((intensities * scaleweights) - (Ih_values * 2.0 * scale_factors * scaleweights))
    dIh = flex.double(np.tile(dIh, n_gbins)) * self.data_manager.g_decay_derivatives
    h_index_cumulative_array = np.array(self.data_manager.h_index_cumulative_array)

    h_index_cumulative_array_tiled = h_index_cumulative_array
    for i in range(1,len(self.data_manager.g_decay)):
      h_index_cumulative_array_tiled = np.append(h_index_cumulative_array_tiled, h_index_cumulative_array[1:]+(i*num))

    dIh = np.add.reduceat(dIh, h_index_cumulative_array_tiled[:-1])/sumgsq_tiled
    dIh = flex.double(np.repeat(dIh, np.tile(self.data_manager.h_index_counter_array, n_gbins)))
    
    drdp = -((Ih_values_tiled * self.data_manager.g_decay_derivatives) + (scale_factors_tiled * dIh))
    grad = (2.0 * rhl_tiled * scaleweights_tiled * drdp) 
    gradient = np.add.reduceat(grad, gradient_reducer)
    return gradient


  def calculate_gradient_fd(self):
    '''calculate gradient array with finite difference approach'''
    delta = 1.0e-6
    gradients = flex.double([0.0]*20)
    for i in range(20):
      intensities = self.data_manager.sorted_reflections['intensity']
      weights = self.data_manager.weights_for_scaling.get_weights()

      for j in range(len(self.data_manager.sorted_reflections['d'])):
        if self.data_manager.sorted_reflections[self.data_manager.active_bin_index][j] == i:
          self.data_manager.sorted_reflections['inverse_scale_factor'][j] -= (0.5*delta)
      self.data_manager.calc_Ih()
      Ih_values_low = self.data_manager.sorted_reflections['Ih_values']
      scale_factors_low = self.data_manager.sorted_reflections['inverse_scale_factor']
      R_low = ((((intensities - (scale_factors_low * Ih_values_low))**2) * weights))

      for j in range(len(self.data_manager.sorted_reflections['d'])):
        if self.data_manager.sorted_reflections[self.data_manager.active_bin_index][j] == i:
          self.data_manager.sorted_reflections['inverse_scale_factor'][j] += (delta)
      self.data_manager.calc_Ih()
      Ih_values_upper = self.data_manager.sorted_reflections['Ih_values']
      scale_factors_upper = self.data_manager.sorted_reflections['inverse_scale_factor']
      R_upper = ((((intensities - (scale_factors_upper * Ih_values_upper))**2) * weights))

      for j in range(len(self.data_manager.sorted_reflections['d'])):
        if self.data_manager.sorted_reflections[self.data_manager.active_bin_index][j] == i:
          self.data_manager.sorted_reflections['inverse_scale_factor'][j] -= (0.5*delta)
      self.data_manager.calc_Ih()

      gradients[i] += (flex.sum(R_upper) - flex.sum(R_low)) / delta
    return gradients

class xds_target_function_log(target_function):
  '''Subclass that takes a data manager object and returns a residual and
  gradient function for a xds-like scaling parameterisation.'''
  def calculate_gradient(self):
    sigma = 0.05
    intensities = self.data_manager.sorted_reflections['intensity']
    variances = self.data_manager.sorted_reflections['variance']
    scale_factors = self.data_manager.sorted_reflections['inverse_scale_factor']
    # list(scale_factors)[100:120]
    Ih_values = self.data_manager.sorted_reflections['Ih_values']
    scaleweights = self.data_manager.weights_for_scaling.get_weights()
    gsq = ((scale_factors)**2) / variances
    sumgsq = np.add.reduceat(gsq, self.data_manager.h_index_cumulative_array[:-1])
    sumgsq = flex.double(np.repeat(sumgsq, self.data_manager.h_index_counter_array))
    dIhdg = ((scale_factors*intensities) - (Ih_values * 2.0 * scale_factors)) / (variances * sumgsq)
    rhl = intensities - (Ih_values * scale_factors)
    grad = ((-2.0 * rhl * ((Ih_values +  dIhdg)*scale_factors) * scaleweights))
    return flex.double(np.bincount(self.data_manager.sorted_reflections[
      self.data_manager.active_bin_index], weights=grad,
      minlength=self.data_manager.active_param_size))

class aimless_target_function(target_function):
  '''Class that takes a data manager object and returns a residual and
  gradient function for an aimless-like scaling parameterisation.'''
  def calculate_gradient(self):
    '''for aimless type scaling, the data manager will need to make use of
    the scale factor derivatives supplied by the basis function'''
    assert 0, 'aimless gradient function not yet defined'
