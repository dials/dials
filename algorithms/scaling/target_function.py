"""
This file defines targets for scaling.

These are initialised with a scaler and an active parameter manager,
and have implementations of residual/gradient calculations for
scaling.
"""
import logging
from dials.array_family import flex
from dials_scaling_helpers_ext import row_multiply
from scitbx import sparse
from dials_scratch_scaling_ext import elementwise_square


class ScalingTarget(object):
  """
  A class to be used by a Scaling Refinery to calculate gradients,
  residuals etc required by the Refinery for minimisation.
  '"""

  _grad_names = ['dI_dp']
  rmsd_names = ["RMSD_I"]
  rmsd_units = ["a.u"]

  def __init__(self, scaler, apm, curvatures=False):
    self.scaler = scaler
    self.apm = apm
    self.weights = self.scaler.Ih_table.weights
    self.curvatures = curvatures

    # Quantities to cache each step
    self._rmsds = None

  def predict(self):#do basis function calcuation and update Ih
    self.scaler.update_for_minimisation(self.apm, self.curvatures)

  def get_num_matches(self):
    """Return the number of reflections in current minimisation cycle."""
    return self.scaler.Ih_table.size

  def rmsds(self):
    """Calculate unweighted RMSDs for the matches."""
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
    """Return the residual vector."""
    Ih_tab = self.scaler.Ih_table
    R = ((((Ih_tab.intensities - (Ih_tab.inverse_scale_factors * Ih_tab.Ih_values))**2)
          * Ih_tab.weights))
    return R

  def calculate_gradients(self):
    """Return a gradient vector on length len(self.apm.x)."""
    Ih_tab = self.scaler.Ih_table
    gsq = ((Ih_tab.inverse_scale_factors)**2) * Ih_tab.weights
    sumgsq = gsq * Ih_tab.h_index_matrix
    rhl = Ih_tab.intensities - (Ih_tab.Ih_values * Ih_tab.inverse_scale_factors)
    dIh = ((Ih_tab.intensities - (Ih_tab.Ih_values * 2.0 *
      Ih_tab.inverse_scale_factors)) * Ih_tab.weights)
    dIh_g = row_multiply(self.apm.derivatives, dIh)
    dIh_red = dIh_g.transpose() * Ih_tab.h_index_matrix
    dIh_by_dpi = row_multiply(dIh_red.transpose(), 1.0/sumgsq)
    term_1 = (-2.0 * rhl * Ih_tab.weights * Ih_tab.Ih_values *
              self.apm.derivatives)
    term_2 = (-2.0 * rhl * Ih_tab.weights * Ih_tab.inverse_scale_factors *
              Ih_tab.h_index_matrix) * dIh_by_dpi
    gradient = term_1 + term_2
    return gradient

  def calculate_curvatures(self):
    """Calculate the second derivatives"""
    # Make some shorthand notation
    Ih_tab = self.scaler.Ih_table
    w = Ih_tab.weights
    I = Ih_tab.intensities
    Ih = Ih_tab.Ih_values
    g = Ih_tab.inverse_scale_factors
    h_idx_transpose = Ih_tab.h_index_matrix.transpose()
    gprime = self.apm.derivatives

    #First calculate Ih' = ((suml wIg') - Ih(suml 2wgg'))/(suml wgg)
    dIh = w * (I - (Ih * 2.0 * g))
    gsq = g * g * w
    v_inv = 1.0/(gsq * Ih_tab.h_index_matrix) #len(n_unique_groups)
    Ihprime_numerator = h_idx_transpose * row_multiply(gprime, dIh)
    Ihprime = row_multiply(Ihprime_numerator, v_inv) #n_unique_groups x n_params matrix

    Ihprime_expanded = (Ihprime.transpose() * Ih_tab.h_expand_matrix).transpose()
    a = row_multiply(gprime, Ih)
    b = row_multiply(Ihprime_expanded, g)
    for j, col in enumerate(a.cols()):
      b[:,j] += col #n_refl x n_params
    rprime = b
    rprimesq = elementwise_square(rprime)
    # First term is sum_refl 2.0 * (r')2
    first_term = 2.0 * w * rprimesq
    print(list(first_term))

    # Second term is sum_refl -2.0 * r * r''
    #  = -2.0 * r * (Ih g'' + 2g'Ih' + g Ih'') = A + B + C
    # Ih '' = u'/v - Ih v'/v. so C = -2 r g u'/v + 2 r g Ih' v'/v = D + E
    r = I - (g * Ih)
    if self.apm.curvatures:
      A = -2.0 * w * r * Ih * self.apm.curvatures #len n_params
    B = row_multiply(gprime, (-4.0 * w * r)).transpose() * Ih_tab.h_index_matrix
    B = B.transpose() * Ihprime_expanded #problem - need to do elementwise multiplication of two big matrices.

    reduced_prefactor = (-2.0 * r * g * Ih_tab.h_index_matrix)
    B = -4.0 * r * ((h_idx_transpose * gprime).transpose() * Ihprime) #len n_params
    vprime_over_v = row_multiply((h_idx_transpose * row_multiply(gprime, (2.0 * w * g))), v_inv) #len n_unique_groups x n_params
    E = (-1.0 * reduced_prefactor * Ihprime) * vprime_over_v

    #E = (2.0 * r * g * Ih_tab.h_index_matrix) * row_multiply(Ihprime, vprime_over_v) #len n_params

    # n_unique_groups x n_params matrices
    if self.apm.curvatures:
      u1 = h_idx_transpose * row_multiply(self.apm.curvatures, w * I)
      u4 = h_idx_transpose * row_multiply(self.apm.curvatures, 2.0 * g * Ih * w)
    u2 = row_multiply((h_idx_transpose * row_multiply(gprime, 2.0 * w * g)), v_inv)
    u3 = h_idx_transpose * row_multiply(elementwise_square(gprime), 2.0 * Ih * w)

    if self.apm.curvatures:
      D1 = reduced_prefactor * row_multiply(u1, v_inv) #len n_params
      D4 = reduced_prefactor * row_multiply(u4, v_inv)
    D2 = (reduced_prefactor * Ihprime) * u2 #len n_params
    D3 = reduced_prefactor * row_multiply(u3, v_inv) #len n_params

    if self.apm.curvatures:
      curvs = A + B + D1 + D2 + D3 + D4 + E + first_term
      print(list(A))
      print(list(B))
      print(list(D1))
      print(list(D2)) #bad
      print(list(D3))
      print(list(D4))
      print(list(E)) #bad
      print(list(first_term))
      print(list(curvs))
      return curvs
    curvs = B + D2 + D3 + E + first_term
    print(list(B))
    print(list(D2))
    print(list(D3))
    print(list(E))
    print(list(first_term))
    print(list(curvs))
    return curvs

  def calculate_jacobian(self):
    """Calculate the jacobian matrix, size Ih_table.size by len(self.apm.x)."""
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
    return flex.sum(self.calculate_residuals()), self.calculate_gradients(), self.calculate_curvatures()

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
  """A special implementation of scaling target for when the scaling is to be
  done against a fixed reference Ih set (i.e scaler is a TargetScaler)
  """

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
