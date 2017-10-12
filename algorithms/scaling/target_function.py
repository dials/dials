'''
Classes that take in a data manager object and return a residual and gradient
vector for minimisation for different scaling parameterisations.
'''
import numpy as np
import copy
from dials.array_family import flex

class target_function(object):
  '''Class that takes in a data manager object and returns a residual
  and gradient function for minimisation.'''
  def __init__(self, data_manager_object):
    self.data_manager = data_manager_object
    #self.parameters = parameters

  def calculate_residual(self):
    '''returns a residual vector'''
    intensities = self.data_manager.reflections_for_scaling['intensity']
    scale_factors = self.data_manager.reflections_for_scaling['inverse_scale_factor']
    Ih_values = self.data_manager.reflections_for_scaling['Ih_values']
    weights = self.data_manager.weights_for_scaling.get_weights()
    R = ((((intensities - (scale_factors * Ih_values))**2) * weights))
    return R

  def calculate_gradient(self):
    '''returns a gradient vector'''
    gradient = flex.double([])
    intensities = self.data_manager.reflections_for_scaling['intensity']
    scale_factors = self.data_manager.reflections_for_scaling['inverse_scale_factor']
    Ih_values = self.data_manager.reflections_for_scaling['Ih_values']
    scaleweights = self.data_manager.weights_for_scaling.get_weights()
    gsq = ((scale_factors)**2) * scaleweights#/ variances
    sumgsq = np.add.reduceat(gsq, self.data_manager.h_index_cumulative_array[:-1])
    rhl = intensities - (Ih_values * scale_factors)
    num = len(intensities)
    dIh = ((intensities * scaleweights) - (Ih_values * 2.0 * scale_factors * scaleweights))

    for i in range(self.data_manager.n_active_params):
      dIh_g = (dIh * self.data_manager.active_derivatives[i*num:(i+1)*num])
      dIh_g = np.add.reduceat(dIh_g, self.data_manager.h_index_cumulative_array[:-1])/sumgsq
      dIh_g = flex.double(np.repeat(dIh_g, self.data_manager.h_index_counter_array))
      drdp = -((Ih_values * self.data_manager.active_derivatives[i*num:(i+1)*num])
               + (scale_factors * dIh_g))
      grad = (2.0 * rhl * scaleweights * drdp)
      gradient.append(flex.sum(grad))
    return gradient

  def return_targets(self):
    '''return residual and gradient arrays'''
    #self.data_manager.update_for_minimisation(self.parameters)
    return self.calculate_residual(), self.calculate_gradient()

  
class xds_target_function_log(target_function):
  '''Subclass that takes a data manager object and returns a residual and
  gradient function for a xds-like scaling parameterisation.'''
  def calculate_gradient(self):
    gradient = flex.double([])
    intensities = self.data_manager.reflections_for_scaling['intensity']
    variances = self.data_manager.reflections_for_scaling['variance']
    scale_factors = self.data_manager.reflections_for_scaling['inverse_scale_factor']
    # list(scale_factors)[100:120]
    Ih_values = self.data_manager.reflections_for_scaling['Ih_values']
    scaleweights = self.data_manager.weights_for_scaling.get_weights()
    gsq = ((scale_factors)**2) *scaleweights
    sumgsq = np.add.reduceat(gsq, self.data_manager.h_index_cumulative_array[:-1])
    #sumgsq = flex.double(np.repeat(sumgsq, self.data_manager.h_index_counter_array))

    rhl = intensities - (Ih_values * scale_factors)
    num = len(intensities)

    dIh = ((scale_factors * intensities) - (Ih_values * 2.0 * scale_factors)) * scaleweights#/ (variances * sumgsq)

    for i in range(self.data_manager.n_active_params):
      dIh_g = (dIh * self.data_manager.active_derivatives[i*num:(i+1)*num])
      dIh_g = np.add.reduceat(dIh_g, self.data_manager.h_index_cumulative_array[:-1])/sumgsq
      dIh_g = flex.double(np.repeat(dIh_g, self.data_manager.h_index_counter_array))
      drdp = -((Ih_values + dIh_g) * scale_factors 
                * self.data_manager.active_derivatives[i*num:(i+1)*num])
      grad = (2.0 * rhl * scaleweights * drdp)
      gradient.append(flex.sum(grad))
    return gradient

