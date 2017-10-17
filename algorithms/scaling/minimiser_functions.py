'''
Classes to create minimiser objects.
'''

import numpy as np
from dials_array_family_flex_ext import *
from cctbx.array_family.flex import *
from cctbx.array_family import flex
from scitbx import lbfgs
from math import exp
import time
#note: include math exp import after flex imports to avoid exp conflicts

class LBFGS_optimiser(object):
  '''Class that takes in Data_Manager object and runs
  an LBFGS minimisation in a Kabsch approach '''
  def __init__(self, Data_Manager_object, param_name):
    self.data_manager = Data_Manager_object
    self.x = None
    self.parameter_name = param_name
    self.data_manager.active_parameterisation = param_name
    self.set_up_parameterisation()
    self.residuals = []
    print "performing scaling on %s reflections out of %s total reflections" % (
      len(self.data_manager.reflections_for_scaling), len(self.data_manager.sorted_reflections))
    self.core_params = lbfgs.core_parameters(maxfev=15)
    self.termination_params = lbfgs.termination_parameters(max_iterations=15)
    lbfgs.run(target_evaluator=self, core_params=self.core_params,
              termination_params=self.termination_params)
    if self.data_manager.scaling_options['parameterization'] == 'standard':
      self.make_all_scales_positive()
    if param_name == 'g_decay':
      if self.data_manager.scaling_options['decay_correction_rescaling']:
        if self.data_manager.scaling_options['parameterization'] == 'standard':
          self.data_manager.scale_gvalues()

  def make_all_scales_positive(self):
    '''catcher that checks all the scale factors are positive in the standard
    parameterization. If they are not, the assumption is that the algorithm
    has got stuck in a local minimum, likely due to a few bad datapoints.
    To cure, the absolute values of the scale factors are taken and the
    minimizer is called again until only positive scale factors are obtained.'''
    if (self.x < 0.0).count(True) > 0.0:
      print """%s of the scale factors is/are negative, taking the absolute 
      values and trying again""" % ((self.x < 0.0).count(True))
      self.x = abs(self.x)
      lbfgs.run(target_evaluator=self, core_params=self.core_params,
                termination_params=self.termination_params)
      if (self.x < 0.0).count(True) > 0.0:
        self.make_all_scales_positive()
      else:
        print "all scales should now be positive"
        self.data_manager.g_parameterisation[self.parameter_name]['parameterisation'] = self.x
    else:
      print "all scales are positive"
      self.data_manager.g_parameterisation[self.parameter_name]['parameterisation'] = self.x

  def return_data_manager(self):
    '''return data_manager method'''
    return self.data_manager

  def set_up_parameterisation(self):
    '''Set up the problem by indicating which g values are being minimised'''
    constant_g_values = []
    for p_type, parameterisation in self.data_manager.g_parameterisation.iteritems():
      bin_index = parameterisation['index']
      if self.parameter_name == p_type:
        self.data_manager.active_bin_index = bin_index
        self.x = parameterisation['parameterisation']
        self.data_manager.n_active_params = len(self.x)
      else:
        constant_g_values.append(self.data_manager.g_parameterisation[p_type]['values'])
    constant_g_values = np.array(constant_g_values)
    if self.data_manager.scaling_options['parameterization'] == 'standard':
      self.data_manager.constant_g_values = flex.double(np.prod(constant_g_values, axis=0))
    elif self.data_manager.scaling_options['parameterization'] == 'log':
      self.data_manager.constant_g_values = flex.double(np.sum(constant_g_values, axis=0))

  def compute_functional_and_gradients(self):
    '''first calculate the updated values of the scale factors and Ih,
    before calculating the residual and gradient functions'''
    self.data_manager.update_for_minimisation(parameters=self.x)
    f, g = self.data_manager.get_target_function()
    f = flex.sum(f)
    self.residuals.append(f)
    print "Residual sum: %12.6g" % f
    return f, g

class aimless_LBFGS_optimiser(object):
  '''Class that takes in Data_Manager object and runs
  an LBFGS minimisation in a Kabsch approach '''
  def __init__(self, Data_Manager_object):
    self.data_manager = Data_Manager_object
    self.x = self.data_manager.active_parameters
    #self.set_up_parameterisation()
    self.residuals = []
    core_params = lbfgs.core_parameters(maxfev=15)
    termination_params = lbfgs.termination_parameters(max_iterations=15)#, traditional_convergence_test_eps=1e2)
    lbfgs.run(target_evaluator=self, core_params=core_params, termination_params = termination_params)

  def return_data_manager(self):
    '''return data_manager method'''
    return self.data_manager

  def compute_functional_and_gradients(self):
    '''first calculate the updated values of the scale factors and Ih,
    before calculating the residual and gradient functions'''
    self.data_manager.update_for_minimisation(parameters=self.x)
    f, g = self.data_manager.get_target_function()
    f = flex.sum(f)
    self.residuals.append(f)
    print "Residual sum: %12.6g" % f
    return f, g


class B_optimiser(object):
  '''Class that takes in minimised decay modulation array and fits a
  global scale and B-factor to the first time row '''
  def __init__(self, Data_Manager_object, initial_values):
    self.data_manager = Data_Manager_object
    d_bin_boundaries = self.data_manager.bin_boundaries['d']
    self.res_values = flex.double([])
    for i in range(0, len(d_bin_boundaries) - 1):
      self.res_values.append(((1.0 / (d_bin_boundaries[i]**2))
                              +(1.0 / (d_bin_boundaries[i+1]**2))) / 2.0)
    self.x = initial_values
    lbfgs.run(target_evaluator=self)

  def compute_functional_and_gradients(self):
    f = self.residual()
    g = self.gradient()
    return f, g

  def residual(self):
    gvalues = self.data_manager.g_decay[0:self.data_manager.binning_parameters['n_d_bins']]
    resolution = self.res_values
    R = 0.0
    for i, val in enumerate(resolution):
      R += ((gvalues[i] * exp(self.x[0] * val)) - self.x[1])**2
    return R

  def gradient(self):
    gvalues = self.data_manager.g_decay[0:self.data_manager.binning_parameters['n_d_bins']]
    resolution = self.res_values
    G = flex.double([0.0, 0.0])
    for i, val in enumerate(resolution):
      G[0] += (2.0 * ((gvalues[i] * exp((self.x[0]) * val)) - self.x[1])
               * resolution[i]*gvalues[i]*exp((self.x[0])*resolution[i]))
      G[1] += -2.0 * ((gvalues[i] * exp((self.x[0]) * val)) - self.x[1])
    return G
