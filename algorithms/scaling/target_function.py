'''
Classes that take in a data manager object and return a residual and gradient
vector for minimisation for different scaling parameterisations.
'''
import numpy as np
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
    sigma = 0.05
    intensities = self.data_manager.sorted_reflections['intensity']
    scale_factors = self.data_manager.sorted_reflections['inverse_scale_factor']
    Ih_values = self.data_manager.sorted_reflections['Ih_values']
    weights = self.data_manager.weights_for_scaling.get_weights()
    #weights = self.data_manager.weights
    '''R = (((((intensities - (scale_factors * Ih_values))**2) / variances)
          + (((scale_factors - 1.0)**2)/(sigma**2))) * weights)'''
    R = ((((intensities - (scale_factors * Ih_values))**2) * weights)
         + (((scale_factors - 1.0)**2)/(sigma**2)))
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
    sigma = 0.05
    intensities = self.data_manager.sorted_reflections['intensity']
    variances = self.data_manager.sorted_reflections['variance']
    scale_factors = self.data_manager.sorted_reflections['inverse_scale_factor']
    Ih_values = self.data_manager.sorted_reflections['Ih_values']
    weights = self.data_manager.weights_for_scaling.get_weights()
    gsq = ((scale_factors)**2) / variances
    sumgsq = np.add.reduceat(gsq, self.data_manager.h_index_cumulative_array[:-1])
    sumgsq = flex.double(np.repeat(sumgsq, self.data_manager.h_index_counter_array))
    dIhdg = flex.double((intensities - (Ih_values * 2.0 * scale_factors))
                        / (variances * sumgsq))
    rhl = flex.double(intensities - (Ih_values * scale_factors))
    grad = ((-2.0 * rhl * ((Ih_values + (scale_factors * dIhdg))) * weights)
            + (2.0 * (scale_factors - 1.0)/(sigma**2)))
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
