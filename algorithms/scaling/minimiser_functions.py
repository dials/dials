'''
Classes to create minimiser objects.
'''
from __future__ import print_function
import numpy as np
from dials_array_family_flex_ext import *
from cctbx.array_family.flex import *
from cctbx.array_family import flex
from scitbx import lbfgs
from math import exp
import time
from data_manager_functions import active_parameter_manager, multi_active_parameter_manager
#note: include math exp import after flex imports to avoid exp conflicts?

class LBFGS_optimiser(object):
  '''Class that takes in Data_Manager object and runs an LBFGS minimisation'''
  def __init__(self, Data_Manager_object, param_name):
    print('\n'+'*'*40+'\n'+'Initialising LBFGS optimiser instance. \n')
    self.data_manager = Data_Manager_object
    if self.data_manager.scaling_options['multi_mode']:
      self.apm = multi_active_parameter_manager(self.data_manager, param_name)
    else:
      self.apm = active_parameter_manager(self.data_manager, param_name)
    self.x = self.apm.x
    print(len(self.x))
    self.residuals = []
    self.core_params = lbfgs.core_parameters(maxfev=15)
    self.termination_params = lbfgs.termination_parameters(max_iterations=15)
    lbfgs.run(target_evaluator=self, core_params=self.core_params,
              termination_params=self.termination_params)
    if param_name == 'g_decay':
      if self.data_manager.scaling_options['decay_correction_rescaling']:
        if self.data_manager.scaling_options['parameterization'] == 'standard':
          self.data_manager.scale_gvalues()
    print(('\nCompleted minimisation for following corrections: {0}\n'
      +'*'*40+'\n').format(''.join(i.lstrip('g_')+' ' for i in param_name)))

  def compute_functional_and_gradients(self):
    '''first calculate the updated values of the scale factors and Ih,
    before calculating the residual and gradient functions'''
    self.data_manager.update_for_minimisation(self.apm)
    f, g = self.data_manager.get_target_function(self.apm)
    f = flex.sum(f)
    self.residuals.append(f)
    print("Residual sum: %12.6g" % f)
    return f, g

  def return_data_manager(self):
    '''return data_manager method'''
    return self.data_manager

  def make_all_scales_positive(self, param_name):
    '''catcher that checks all the scale factors are positive in the standard
    parameterization. If they are not, the assumption is that the algorithm
    has got stuck in a local minimum, likely due to a few bad datapoints.
    To cure, the absolute values of the scale factors are taken and the
    minimizer is called again until only positive scale factors are obtained.'''
    if (self.x < 0.0).count(True) > 0.0:
      print("""%s of the scale factors is/are negative, taking the absolute
      values and trying again""" % ((self.x < 0.0).count(True)))
      self.x = abs(self.x)
      lbfgs.run(target_evaluator=self, core_params=self.core_params,
                termination_params=self.termination_params)
      if (self.x < 0.0).count(True) > 0.0:
        self.make_all_scales_positive(param_name)
      else:
        print("all scales should now be positive")
    else:
      print("all scales are positive")

class error_scale_LBFGSoptimiser(object):
  def __init__(self, Ih_table, starting_values):
    # default start a = 1.0, b = 0.0, c = 0.0
    self.Ih_table = Ih_table
    self.x = starting_values
    self.Ih_table.Ih_table['sigmaprime'] = self.calc_sigmaprime()
    self.Ih_table.Ih_table['delta_hl'] = self.calc_deltahl()
    self.bin_intensities()

    lbfgs.run(target_evaluator=self)
    print("minimised error scales, values are %s" % list(self.x))
    #import matplotlib.pyplot as plt
    #plt.hist(self.Ih_table.Ih_table['delta_hl'], 60)
    #plt.show()
    #plt.plot(self.bin_vars)
    #plt.show()

  def compute_functional_and_gradients(self):
    '''first calculate the updated values of sigmaprime and delta_hl,
    before calculating the residual and gradient functions'''
    self.Ih_table.Ih_table['sigmaprime'] = self.calc_sigmaprime()
    self.Ih_table.Ih_table['delta_hl'] = self.calc_deltahl()
    R = self.calc_residual()
    G = self.calc_gradient()
    return R, G

  def calc_sigmaprime(self):
    sigmaprime = self.x[0] * ((1.0/self.Ih_table.Ih_table['weights'])
      #+ (self.x[1]*self.Ih_table.Ih_table['intensity'])
      + ((self.x[1]*self.Ih_table.Ih_table['intensity'])**2))**0.5
    return sigmaprime

  def calc_deltahl(self):
    #first calcualte n from h_index_cumulative_array?
    n_h = self.Ih_table.n_h
    I_hl = self.Ih_table.Ih_table['intensity']
    g_hl = self.Ih_table.Ih_table['inverse_scale_factor']
    I_h = self.Ih_table.Ih_table['Ih_values']
    prefactor = ((n_h - flex.double([1.0]*len(n_h))) / n_h)**0.5
    delta_hl = prefactor * ((I_hl/g_hl) - I_h) / self.Ih_table.Ih_table['sigmaprime']
    return delta_hl

  def bin_intensities(self):
    sel = flex.sort_permutation(self.Ih_table.Ih_table['intensity'])
    self.Ih_table.Ih_table = self.Ih_table.Ih_table.select(sel)
    deltahl = self.Ih_table.Ih_table['delta_hl']
    n = len(self.Ih_table.Ih_table)
    n_bins = 20
    self.n_bin_count = []
    for i in range(0, n_bins+1):
      self.n_bin_count.append((i*n)//n_bins)
    self.n_in_each_bin = flex.double([])
    for i, val in enumerate(self.n_bin_count[:-1]):
      self.n_in_each_bin.append(self.n_bin_count[i+1]-val)

  def calc_residual(self):
    deltahl = self.Ih_table.Ih_table['delta_hl']
    import numpy as np
    sum_deltasq = flex.double(np.add.reduceat(deltahl**2, self.n_bin_count[:-1]))
    sum_delta_sq = flex.double(np.add.reduceat(deltahl, self.n_bin_count[:-1]))**2
    self.bin_vars = ((sum_deltasq/flex.double(self.n_in_each_bin)) -
                     (sum_delta_sq/(flex.double(self.n_in_each_bin)**2)))
    #print list(self.bin_vars)
    R = flex.sum(((flex.double([1.0]*len(self.bin_vars)) - self.bin_vars)**2))
    R = R + 1.0*((1.0 - self.x[0])**2 + (self.x[1]**2))# + (self.x[2]**2))
    return R

  def calc_gradient(self):
    I_hl = self.Ih_table.Ih_table['intensity']
    sigmaprime = self.Ih_table.Ih_table['sigmaprime']
    delta_hl = self.Ih_table.Ih_table['delta_hl']
    dsig_da = sigmaprime/self.x[0]
    #dsig_db = I_hl * (self.x[0]**2) / (2.0 * sigmaprime)
    dsig_dc = self.x[1] * (I_hl**2) * (self.x[0]**2) / sigmaprime
    ddelta_dsigma = -1.0 * delta_hl / sigmaprime
    dsig_list = [ddelta_dsigma * dsig_da,# ddelta_dsigma * dsig_db,
                 ddelta_dsigma * dsig_dc]
    gradient = flex.double([])
    for deriv in dsig_list:
      term1 = flex.double(np.add.reduceat(2.0 * delta_hl * deriv, self.n_bin_count[:-1]))
      term2a = flex.double(np.add.reduceat(delta_hl, self.n_bin_count[:-1]))
      term2b = flex.double(np.add.reduceat(deriv, self.n_bin_count[:-1]))
      grad = -1.0 * (2.0 * (flex.double([1.0]*len(self.bin_vars)) - self.bin_vars) *
              ((term1 / self.n_in_each_bin)
               - (2.0 * term2a * term2b / (self.n_in_each_bin**2))))
      gradient.append(flex.sum(grad))
    gradient = gradient + flex.double([-1.0 * 2.0 * (1.0 - self.x[0]),
                                       1.0 * 2.0 * self.x[1]])
    return gradient


