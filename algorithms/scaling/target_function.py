'''
Classes that take in a scaler and return a residual and gradient
vector for minimisation for different scaling parameterisations.
'''
import abc
from dials.array_family import flex
from dials_scaling_helpers_ext import row_multiply
from scitbx import sparse


class ScalingTarget(object):
  '''A class to be used by a Scaling Refinery to calculate gradients,
  residuals etc required by the Refinery for minimisation.
  '''
  __metaclass__ = abc.ABCMeta
  _grad_names = ['dI_dp']
  rmsd_names = ["RMSD_I"]
  rmsd_units = ["a.u"]

  def __init__(self, scaler, apm):
    self.scaler = scaler
    self.apm = apm
    self.weights = self.scaler.Ih_table.weights

    # Quantities to cache each step
    self._rmsds = None

  def predict(self):#do basis function calcuation and update Ih
    self.scaler.update_for_minimisation(self.apm)

  def get_num_matches(self):
    return self.scaler.Ih_table.size

  def rmsds(self):
    """calculate unweighted RMSDs for the matches"""
    # cache rmsd calculation for achieved test
    R = self.calculate_residuals()
    if 'absorption' in self.apm.components_list:
      restr = self.scaler.calc_absorption_restraint(self.apm)
      if restr:
        R.extend(restr[0])
    self._rmsds = [(flex.sum((R))/self.scaler.Ih_table.size)**0.5]
    #print("rmsds %s" % self._rmsds)
    return self._rmsds

  def achieved(self):
    return False #implement a method here?

  def calculate_residuals(self):
    '''returns a residual vector'''
    Ih_tab = self.scaler.Ih_table
    R = ((((Ih_tab.intensities - (Ih_tab.inverse_scale_factors * Ih_tab.Ih_values))**2)
          * Ih_tab.weights))
    return R

  def calculate_gradients(self):
    '''returns a gradient vector on length len(self.apm.x)'''
    Ih_tab = self.scaler.Ih_table
    gsq = ((Ih_tab.inverse_scale_factors)**2) * Ih_tab.weights
    sumgsq = gsq * Ih_tab.h_index_matrix
    rhl = Ih_tab.intensities - (Ih_tab.Ih_values * Ih_tab.inverse_scale_factors)
    dIh = ((Ih_tab.intensities - (Ih_tab.Ih_values * 2.0 * Ih_tab.inverse_scale_factors))
           * Ih_tab.weights)
    dIh_g = row_multiply(self.apm.derivatives, dIh)
    dIh_red = dIh_g.transpose() * Ih_tab.h_index_matrix
    dIh_by_dpi = row_multiply(dIh_red.transpose(), 1.0/sumgsq)
    term_1 = (-2.0 * rhl * Ih_tab.weights * Ih_tab.Ih_values *
              self.apm.derivatives)
    term_2 = (-2.0 * rhl * Ih_tab.weights * Ih_tab.inverse_scale_factors *
              Ih_tab.h_index_matrix) * dIh_by_dpi
    gradient = term_1 + term_2
    return gradient

  def calculate_jacobian(self):
    '''returns a jacobian matrix of size Ih_table.size by len(self.apm.x)'''
    Ih_tab = self.scaler.Ih_table
    invsumgsq = 1.0 / ((Ih_tab.inverse_scale_factors)**2) * Ih_tab.h_index_matrix
    dIh = ((Ih_tab.intensities - (Ih_tab.Ih_values * 2.0 * Ih_tab.inverse_scale_factors)))
    dIh_g = row_multiply(self.apm.derivatives, dIh)
    dIh_red = dIh_g.transpose() * Ih_tab.h_index_matrix
    dIh_by_dpi = row_multiply(dIh_red.transpose(), invsumgsq)
    dIh_by_dpi = dIh_by_dpi.transpose() * Ih_tab.h_expand_matrix
    term1 = row_multiply(dIh_by_dpi.transpose(), Ih_tab.inverse_scale_factors)
    term2 = row_multiply(self.apm.derivatives, Ih_tab.Ih_values)
    for i, col in enumerate(term2.cols()):
      term1[:, i] += col # Sum columns to create the jacobian
    return term1 # Return the jacobian

  # The following methods are for adaptlbfgs.
  def compute_functional_gradients(self):
    return flex.sum(self.calculate_residuals()), self.calculate_gradients()

  def compute_functional_gradients_and_curvatures(self):
    return flex.sum(self.calculate_residuals()), self.calculate_gradients(), None #curvatures

  def compute_restraints_functional_gradients_and_curvatures(self):
    'calculate restraints on parameters'
    restraints = None
    if 'absorption' in self.apm.components_list:
      restr = self.scaler.calc_absorption_restraint(self.apm)
      if restr:
        resid_restr = flex.sum(restr[0]) #want just a value to add to total functional here
        grad_restr = restr[1]
        restraints = [resid_restr, grad_restr, None]
    return restraints #list of restraints to add to resid, grads and curvs?

  'the following methods are for adaptlstbx (GN/ LM algorithms)'
  def compute_residuals(self):
    return self.calculate_residuals(), self.weights

  def compute_residuals_and_gradients(self):
    return self.calculate_residuals(), self.calculate_jacobian(), self.weights

  def compute_restraints_residuals_and_gradients(self):
    'calculate restraints on parameters'
    restraints = None
    if 'absorption' in self.apm.components_list:
      restraints = self.scaler.compute_restraints_residuals_jacobian(self.apm)
    return restraints

class ScalingTargetFixedIH(ScalingTarget):
  '''A special implementation of scaling target for when the scaling is to be
  done against a fixed reference Ih set (i.e scaler is a TargetScaler)
  '''
  def calculate_gradients(self):
    G = flex.double([])
    for i, unscaled_scaler in enumerate(self.scaler.unscaled_scalers):
      apm = self.apm.apm_list[i]
      Ih_tab = unscaled_scaler.Ih_table
      rhl = Ih_tab.intensities - (Ih_tab.Ih_values * Ih_tab.inverse_scale_factors)
      G.extend(-2.0 * rhl * Ih_tab.weights * Ih_tab.Ih_values * apm.derivatives)
    return G

  def calculate_residuals(self):
    '''returns a residual vector'''
    R = flex.double([])
    for unscaled_scaler in self.scaler.unscaled_scalers:
      Ih_tab = unscaled_scaler.Ih_table
      R.extend((((Ih_tab.intensities - (Ih_tab.inverse_scale_factors * Ih_tab.Ih_values))**2)
                * Ih_tab.weights))
    return R

  def compute_residuals_and_gradients(self):
    assert 0, 'method not yet implemented for targeted scaling'
