'''
Classes that take in a data manager object and minimisation parameters and
return the scale factors and derivatives of the scale factors w.r.t.
the parameters
'''
from cctbx.array_family import flex
import numpy as np

class basis_function(object):
  '''Superclass for basis function that takes in a data manager object and
  minimisation parameters and calcuates the scale factors and derivatives'''
  def __init__(self, data_manager_object, parameters):
    self.data_manager = data_manager_object
    self.parameters = parameters

  def calculate_scale_factors(self):
    '''Calculate scale factors method to be overloaded by subclass'''
    assert 0, 'parameterization-specific scale factors must be defined in a subclass'

  def calculate_derivatives(self):
    '''Calculate derivatives method to be overloaded by subclass'''
    assert 0, 'parameterization-specific derivatives must be defined in a subclass'

  def return_basis(self):
    '''Return the calculated scale factors and derivatives to be used
    in minimisation'''
    return self.calculate_scale_factors(), self.calculate_derivatives()

class xds_basis_function(basis_function):
  '''Subclass of basis_function for xds parameterisation'''
  def calculate_scale_factors(self):
    gxvalues = flex.double([self.parameters[i] for i in
           self.data_manager.sorted_reflections[self.data_manager.active_bin_index]])
    scale_factors = gxvalues * self.data_manager.constant_g_values
    return scale_factors

  def calculate_derivatives(self):
    '''xds target function does not require derivatives, so None returned'''
    return None

class xds_basis_function_log(basis_function):
  '''Subclass of basis_function for xds parameterisation'''
  def calculate_scale_factors(self):
    gxvalues = flex.double([self.parameters[i] for i in
           self.data_manager.sorted_reflections[self.data_manager.active_bin_index]])
    scale_factors = gxvalues + self.data_manager.constant_g_values
    #print list(scale_factors)[100:120]
    return flex.double(np.exp(scale_factors))

  def calculate_derivatives(self):
    '''xds target function does not require derivatives, so None returned'''
    return None