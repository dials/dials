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
  def __init__(self, data_manager_object, parameters=None):
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
           self.data_manager.reflections_for_scaling[self.data_manager.active_bin_index]])
    scale_factors = gxvalues * self.data_manager.constant_g_values
    return scale_factors

  def calculate_derivatives(self):
    '''Derivatives are fixed by the parameterisation'''
    derivatives = self.data_manager.g_parameterisation[
      self.data_manager.active_parameterisation]['derivatives']
    return derivatives

class aimless_basis_function(basis_function):
  '''Subclass of basis_function for xds parameterisation'''
  def calculate_scale_factors(self):
    ngscale = self.data_manager.n_g_scale_params
    ngdecay = self.data_manager.n_g_decay_params
    d = self.data_manager.reflections_for_scaling['d']
    B = self.data_manager.active_parameters[ngscale:ngscale+ngdecay]
    scale = self.data_manager.active_parameters[0:ngscale]
    t2 = self.data_manager.reflections_for_scaling['t2_bin_index']
    t1 = self.data_manager.reflections_for_scaling['t1_bin_index']
    B_factor = flex.double([np.exp(B[t]/(2*(d[i]**2))) for i,t in enumerate(t2)])
    scale_factor = flex.double([scale[i] for i in t1])
    #for absorption, will want to add a calculation of the scattering vector to the refl table
    absorption_factor = flex.double([1.0 for _ in d])
    inverse_scale_factors = B_factor * scale_factor * absorption_factor
    self.data_manager.reflections_for_scaling['B_factor'] = B_factor
    self.data_manager.reflections_for_scaling['scale_factor'] = scale_factor
    self.data_manager.reflections_for_scaling['absorption_factor'] = absorption_factor
    return inverse_scale_factors

  def calculate_derivatives(self):
    '''xds target function does not require derivatives, so None returned'''
    total_derivatives = []
    ngscale = self.data_manager.n_g_scale_params
    ngdecay = self.data_manager.n_g_decay_params
    ngabs = self.data_manager.n_g_abs_params
    n = len(self.data_manager.reflections_for_scaling)
    gscale_derivatives = flex.double([0.0]* ngscale * n)
    for i, l in enumerate(self.data_manager.reflections_for_scaling['t1_bin_index']):
      gscale_derivatives[(l*n)+i] = 1.0
    total_derivatives += (flex.double(np.repeat((self.data_manager.reflections_for_scaling['B_factor'] 
                                                 * self.data_manager.reflections_for_scaling['absorption_factor']), 
                                                ngscale)) * gscale_derivatives)
    gdecay_derivatives = flex.double([0.0] * ngdecay * n)
    d = self.data_manager.reflections_for_scaling['d']
    B = self.data_manager.active_parameters[ngscale:ngscale+ngdecay]
    for i, l in enumerate(self.data_manager.reflections_for_scaling['t2_bin_index']):
      gdecay_derivatives[(l*n)+i] = np.exp(B[l] / (2 * (d[i]**2))) / (2 * (d[i]**2))
    total_derivatives += (flex.double(np.repeat((self.data_manager.reflections_for_scaling['scale_factor'] 
                                                 * self.data_manager.reflections_for_scaling['absorption_factor']), 
                                                ngdecay)) * gdecay_derivatives)
    gabs_derivatives = flex.double([0.0]* ngabs * n)
    total_derivatives += (flex.double(np.repeat((self.data_manager.reflections_for_scaling['scale_factor'] 
                                                 * self.data_manager.reflections_for_scaling['B_factor']), 
                                                ngabs)) * gabs_derivatives)
    return flex.double(total_derivatives)

class xds_basis_function_log(basis_function):
  '''Subclass of basis_function for xds parameterisation'''
  def calculate_scale_factors(self):
    gxvalues = flex.double([self.parameters[i] for i in
           self.data_manager.reflections_for_scaling[self.data_manager.active_bin_index]])
    scale_factors = gxvalues + self.data_manager.constant_g_values
    #print list(scale_factors)[100:120]
    return flex.double(np.exp(scale_factors))

  def calculate_derivatives(self):
    '''Derivatives are fixed by the parameterisation'''
    derivatives = self.data_manager.g_parameterisation[
      self.data_manager.active_parameterisation]['derivatives']
    return derivatives
