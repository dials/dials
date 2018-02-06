'''
Classes to create minimiser objects.
'''
from __future__ import print_function
import logging
from dials.array_family import flex
from scitbx import lbfgs, sparse
from dials.algorithms.scaling.ParameterHandler import \
  multi_active_parameter_manager, active_parameter_manager, target_multi_active_parameter_manager

logger = logging.getLogger('dials')

class LBFGS_optimiser(object):
  '''Class that takes in scaler object and runs an LBFGS minimisation'''
  def __init__(self, scaler, param_lists):
    logger.info(('\n'+'*'*40+'\n'+'Initialising LBFGS optimiser instance. \n'))
    self.scaler = scaler
    from dials.algorithms.scaling.ScalerFactory import MultiScaler, TargetScaler
    if isinstance(self.scaler, TargetScaler) and len(self.scaler.unscaled_scalers) > 1:
      self.apm = target_multi_active_parameter_manager(self.scaler, param_lists)
    elif isinstance(self.scaler, TargetScaler) and len(self.scaler.unscaled_scalers) == 1:
      self.apm = active_parameter_manager(self.scaler.unscaled_scalers[0], param_lists[0])
    elif isinstance(self.scaler, MultiScaler):
      self.apm = multi_active_parameter_manager(self.scaler, param_lists)
    else:
      self.apm = active_parameter_manager(self.scaler, param_lists[0])
    self.x = self.apm.x
    self.residuals = []
    self.core_params = lbfgs.core_parameters(maxfev=15)
    self.termination_params = lbfgs.termination_parameters(max_iterations=15)
    lbfgs.run(target_evaluator=self, core_params=self.core_params,
              termination_params=self.termination_params)
    #if param_name == 'g_decay':
    #  if self.scaler.params.scaling_options.decay_correction_rescaling:
    #    if self.scaler.params.scaling_options.minimisation_parameterisation == 'standard':
    #      self.scaler.scale_gvalues()
    #logger.info(('\nCompleted minimisation for following corrections: {0}\n'
    #       +'*'*40+'\n').format(''.join(i.lstrip('g_')+' ' for i in param_lists)))

  def compute_functional_and_gradients(self):
    '''first calculate the updated values of the scale factors and Ih,
    before calculating the residual and gradient functions'''
    self.scaler.update_for_minimisation(self.apm)
    f, g = self.scaler.get_target_function(self.apm)
    logger.debug('\nParameter values \n')
    logger.debug(str(list(self.x)) + '\n')
    logger.debug('Parameter derivatives \n')
    logger.debug(str(list(g)) + '\n')
    #f = flex.sum(f)
    self.residuals.append(f)
    logger.info("Residual sum: %12.6g" % f)
    return f, g

  def return_scaler(self):
    '''return scaler method'''
    from dials.algorithms.scaling.ScalerFactory import MultiScaler
    if not isinstance(self.scaler, MultiScaler):
      if 'g_scale' in self.apm.active_parameterisation:
        self.scaler.normalise_scale_component()
      if 'g_decay' in self.apm.active_parameterisation:
        self.scaler.normalise_decay_component()
    return self.scaler

  def make_all_scales_positive(self, param_name):
    '''catcher that checks all the scale factors are positive in the standard
    parameterization. If they are not, the assumption is that the algorithm
    has got stuck in a local minimum, likely due to a few bad datapoints.
    To cure, the absolute values of the scale factors are taken and the
    minimizer is called again until only positive scale factors are obtained.'''
    if (self.x < 0.0).count(True) > 0.0:
      logger.info("""%s of the scale factors is/are negative, taking the absolute
      values and trying again""" % ((self.x < 0.0).count(True)))
      self.x = abs(self.x)
      lbfgs.run(target_evaluator=self, core_params=self.core_params,
                termination_params=self.termination_params)
      if (self.x < 0.0).count(True) > 0.0:
        self.make_all_scales_positive(param_name)
      else:
        logger.info("all scales should now be positive")
    else:
      logger.info("all scales are positive")

class error_scale_LBFGSoptimiser(object):
  '''Class that minimises an error model for an Ih_table'''
  def __init__(self, Ih_table, starting_values):
    # default start a = 1.0, b = 0.05
    # note - don't initialise with b(SdAdd) = 0.0 or it gets stuck on 0!!
    self.Ih_table = Ih_table
    self.x = starting_values
    self.sigmaprime = None
    self.delta_hl = None
    self.calc_sigmaprime()
    self.calc_deltahl()
    self.bin_intensities()
    self.bin_vars = None
    logger.info("Initialised error model LBFGS optimiser instance. \n")
    lbfgs.run(target_evaluator=self)
    logger.info("Minimised error model with parameters {0:.5f} and {1:.5f}. {sep}"
          .format(self.x[0], self.x[1], sep='\n'))

  def compute_functional_and_gradients(self):
    '''first calculate the updated values of sigmaprime and delta_hl,
    before calculating the residual and gradient functions'''
    self.calc_sigmaprime()
    self.calc_deltahl()
    R = self.calc_error_residual()
    G = self.calc_error_gradient()
    return R, G

  def calc_sigmaprime(self):
    '''function to calculate the updated standard deviation'''
    sigmaprime = self.x[0] * ((1.0/self.Ih_table.weights)
      #+ (self.x[1]*self.Ih_table.Ih_table['intensity'])
      + ((self.x[1]*self.Ih_table.intensities)**2))**0.5
    self.sigmaprime = sigmaprime

  def calc_deltahl(self):
    '''function to calculate the normalised deviation of the intensities'''
    n_h = self.Ih_table.n_h
    I_hl = self.Ih_table.intensities
    g_hl = self.Ih_table.inverse_scale_factors
    I_h = self.Ih_table.Ih_values
    prefactor = ((n_h - flex.double([1.0]*len(n_h))) / n_h)**0.5
    delta_hl = prefactor * ((I_hl/g_hl) - I_h) / self.sigmaprime
    self.delta_hl = delta_hl

  def bin_intensities(self):
    '''bin data into intensity bins, and create a 'bin_reducer' matrix for
       summation over indices '''
    sel = flex.sort_permutation(self.Ih_table.intensities)
    #self.Ih_table = self.Ih_table.select(sel)
    self.sigmaprime = self.sigmaprime.select(sel)
    self.delta_hl = self.delta_hl.select(sel)
    n = len(self.delta_hl)
    if n < 10000: # what is a sensible limit here?
      n_bins = 10 # what is a sensible number of bins?
    else:
      n_bins = 20
    self.n_bin_cumulative_array = []
    for i in range(0, n_bins+1):
      self.n_bin_cumulative_array.append((i*n)//n_bins)
    self.n_bin_counter_array = flex.double([])
    for i, val in enumerate(self.n_bin_cumulative_array[:-1]):
      self.n_bin_counter_array.append((self.n_bin_cumulative_array[i+1]-val))
    n = self.n_bin_cumulative_array[-1]
    self.bin_reducer = sparse.matrix(n, len(self.n_bin_counter_array))
    for i in range(len(self.n_bin_cumulative_array)-1):
      col = sparse.matrix_column(n)
      start_idx = self.n_bin_cumulative_array[i]
      for j in range(int(self.n_bin_counter_array[i])):
        col[start_idx+j] = 1
      self.bin_reducer[:, i] = col

  def calc_error_residual(self):
    'calculate the residual'
    deltahl = self.delta_hl
    sum_deltasq = (deltahl**2) * self.bin_reducer
    sum_delta_sq = (deltahl * self.bin_reducer)**2
    self.bin_vars = ((sum_deltasq/flex.double(self.n_bin_counter_array)) -
                     (sum_delta_sq/(flex.double(self.n_bin_counter_array)**2)))
    R = flex.sum(((flex.double([1.0]*len(self.bin_vars)) - self.bin_vars)**2))
    R = R + (25.0*((1.0 - self.x[0])**2)) + (400.0*((0.0001 - self.x[1])**2))
    return R

  def calc_error_gradient(self):
    'calculate the gradient vector'
    I_hl = self.Ih_table.intensities
    sigmaprime = self.sigmaprime
    delta_hl = self.delta_hl
    dsig_da = sigmaprime/self.x[0]
    #dsig_db = I_hl * (self.x[0]**2) / (2.0 * sigmaprime)
    dsig_dc = self.x[1] * (I_hl**2) * (self.x[0]**2) / sigmaprime
    ddelta_dsigma = -1.0 * delta_hl / sigmaprime
    dsig_list = [ddelta_dsigma * dsig_da,# ddelta_dsigma * dsig_db,
                 ddelta_dsigma * dsig_dc]
    gradient = flex.double([])
    for deriv in dsig_list:
      term1 = 2.0 * delta_hl * deriv * self.bin_reducer
      term2a = delta_hl * self.bin_reducer
      term2b = deriv * self.bin_reducer
      grad = (2.0 * (flex.double([1.0]*len(self.bin_vars)) - self.bin_vars)
              * ((term1 / self.n_bin_counter_array)
                 - (2.0 * term2a * term2b / (self.n_bin_counter_array**2)))) * -1.0
      gradient.append(flex.sum(grad))
    gradient = gradient + flex.double([-1.0 * 2.0 * (1.0 - self.x[0]),
                                       1.0 * 2.0 * self.x[1]])
    return gradient
