'''
Classes that take in a data manager object and return a residual and gradient
vector for minimisation for different scaling parameterisations.
'''
import numpy as np
import copy
from dials.array_family import flex
from dials_scaling_helpers_ext import row_multiply
import time as time
from dials_scaling_helpers_ext import row_multiply
from scitbx import sparse

class target_function(object):
  '''Class that takes in a data manager object and returns a residual
  and gradient function for minimisation.'''
  def __init__(self, data_manager_object, apm):
    self.data_manager = data_manager_object
    self.apm = apm

  def calculate_residual(self):
    '''returns a residual vector'''
    intensities = self.data_manager.Ih_table.Ih_table['intensity']
    scale_factors = self.data_manager.Ih_table.Ih_table['inverse_scale_factor']
    Ih_values = self.data_manager.Ih_table.Ih_table['Ih_values']
    weights = self.data_manager.Ih_table.Ih_table['weights']
    R = ((((intensities - (scale_factors * Ih_values))**2) * weights))
    if self.data_manager.scaling_options['scaling_method'] == 'aimless':
      if 'g_absorption' in self.apm.active_parameterisation:
        constraint_values = self.data_manager.calc_absorption_constraint(self.apm)[0]
        R.extend(constraint_values)
    return R

  def calculate_gradient(self):
    '''returns a gradient vector'''
    gradient = flex.double([])
    intensities = self.data_manager.Ih_table.Ih_table['intensity']
    scale_factors = self.data_manager.Ih_table.Ih_table['inverse_scale_factor']
    Ih_values = self.data_manager.Ih_table.Ih_table['Ih_values']
    scaleweights = self.data_manager.Ih_table.Ih_table['weights']
    gsq = ((scale_factors)**2) * scaleweights#/ variances
    sumgsq = flex.double(np.add.reduceat(gsq, self.data_manager.Ih_table.h_index_cumulative_array[:-1]))
    rhl = intensities - (Ih_values * scale_factors)
    num = len(intensities)
    dIh = ((intensities * scaleweights) - (Ih_values * 2.0 * scale_factors * scaleweights))
    
    #dIh = flex.double(np.repeat(dIh, self.apm.n_active_params))
    '''reshaped_derivs = np.reshape(np.array(self.data_manager.active_derivatives), (self.apm.n_active_params, num))
    dIh_g = reshaped_derivs * dIh
    #dIh_g = np.reshape(dIh_g, (self.apm.n_active_params, num))
    dIh_g = dIh_g.transpose()
    dIh_g = np.add.reduceat(dIh_g, self.data_manager.Ih_table.h_index_cumulative_array[:-1])
    dIh_g = dIh_g.transpose()
    dIh_g = dIh_g * 1.0/sumgsq
    dIh_g = np.repeat(dIh_g, np.array(self.data_manager.Ih_table.h_index_counter_array),axis=1)
    drdp = -((reshaped_derivs * Ih_values) + (dIh_g * scale_factors))
    grad = drdp * (2.0 * rhl * scaleweights)
    gradient = flex.double(np.sum(grad, axis=1))'''


    dIh_g = row_multiply(self.data_manager.active_derivatives, dIh)
    dIh_g = dIh_g.transpose()
    dIh_g = dIh_g * self.data_manager.Ih_table.h_index_mat
    dIh_g = dIh_g.transpose()
    #print len(sumgsq)
    #print dIh_g.n_rows
    #print dIh_g.n_cols
    dIh_g = row_multiply(dIh_g, 1.0/sumgsq)
    dIh_by_dpi = dIh_g.transpose()
    
    #n = len(rhl)
    #red_2gr = sparse.matrix(n,1)
    reduced_2gr = 2.0 * rhl * scaleweights * scale_factors
    part_2 = reduced_2gr * self.data_manager.Ih_table.h_index_mat
    part_2 = part_2 * dIh_by_dpi.transpose()

    reduced_2rIh = 2.0 * rhl * scaleweights * Ih_values
    part_1 = reduced_2rIh * self.data_manager.active_derivatives
    total = part_1 + part_2
    print list(total)
    return total
    '''exit()
    

    for i in range(self.apm.n_active_params):
      dIh_g = (dIh * self.data_manager.active_derivatives[i*num:(i+1)*num])
      dIh_g = np.add.reduceat(dIh_g, self.data_manager.Ih_table.h_index_cumulative_array[:-1])/sumgsq
      dIh_g = flex.double(np.repeat(dIh_g, self.data_manager.Ih_table.h_index_counter_array))
      drdp = -((Ih_values * self.data_manager.active_derivatives[i*num:(i+1)*num])
               + (scale_factors * dIh_g))
      grad = (2.0 * rhl * scaleweights * drdp)
      gradient.append(flex.sum(grad))
    if self.data_manager.scaling_options['scaling_method'] == 'aimless':
      if 'g_absorption' in self.apm.active_parameterisation:
        gradient += self.data_manager.calc_absorption_constraint(self.apm)[1]
    return gradient'''

  def return_targets(self):
    '''return residual and gradient arrays'''
    return self.calculate_residual(), self.calculate_gradient()

class target_function_fixedIh(target_function):
  def __init__(self, data_manager_object, apm):
    target_function.__init__(self, data_manager_object, apm)

  'subclass to calculate the gradient for KB scaling against a fixed Ih'
  def calculate_gradient(self):
    gradient = flex.double([])
    intensities = self.data_manager.Ih_table.Ih_table['intensity']
    scale_factors = self.data_manager.Ih_table.Ih_table['inverse_scale_factor']
    Ih_values = self.data_manager.Ih_table.Ih_table['Ih_values']
    scaleweights = self.data_manager.Ih_table.Ih_table['weights']
    rhl = intensities - (Ih_values * scale_factors)
    for i in range(self.apm.n_active_params):
      drdp = -Ih_values * self.data_manager.active_derivatives[i*len(intensities):
                                                               (i+1)*len(intensities)]
      grad = (2.0 * rhl * scaleweights * drdp)
      gradient.append(flex.sum(grad))
    return gradient

class xds_target_function_log(target_function):
  '''Subclass that takes a data manager object and returns a residual and
  gradient function for a xds-like scaling parameterisation.'''
  def calculate_gradient(self):
    gradient = flex.double([])
    intensities = self.data_manager.Ih_table.Ih_table['intensity']
    scale_factors = self.data_manager.Ih_table.Ih_table['inverse_scale_factor']
    Ih_values = self.data_manager.Ih_table.Ih_table['Ih_values']
    scaleweights = self.data_manager.Ih_table.Ih_table['weights']
    gsq = ((scale_factors)**2) *scaleweights
    sumgsq = np.add.reduceat(gsq, self.data_manager.Ih_table.h_index_cumulative_array[:-1])
    #sumgsq = flex.double(np.repeat(sumgsq, self.data_manager.h_index_counter_array))

    rhl = intensities - (Ih_values * scale_factors)
    num = len(intensities)
    dIh = ((scale_factors * intensities) - (Ih_values * 2.0 * scale_factors)) * scaleweights
    for i in range(self.data_manager.n_active_params):
      dIh_g = (dIh * self.data_manager.active_derivatives[i*num:(i+1)*num])
      dIh_g = np.add.reduceat(dIh_g, self.data_manager.Ih_table.h_index_cumulative_array[:-1])/sumgsq
      dIh_g = flex.double(np.repeat(dIh_g, self.data_manager.Ih_table.h_index_counter_array))
      drdp = -((Ih_values + dIh_g) * scale_factors
                * self.data_manager.active_derivatives[i*num:(i+1)*num])
      grad = (2.0 * rhl * scaleweights * drdp)
      gradient.append(flex.sum(grad))
    return gradient
