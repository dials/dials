"""
This file defines targets for scaling.

These are initialised with a scaler and an active parameter manager,
and have implementations of residual/gradient calculations for
scaling.
"""
from dials.array_family import flex
from scitbx import sparse
from dials_scaling_helpers_ext import row_multiply
from dials_scratch_scaling_ext import elementwise_square
from basis_functions import basis_function

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
    if self.scaler.params.scaling_options.use_free_set:
      self.rmsd_names = ["RMSD_I", "Free RMSD_I"]
      self.rmsd_units = ["a.u", "a.u"]
    # Quantities to cache each step
    self._rmsds = None

  def predict(self):
    """Calculate the current values of quantities for the minimisation. Defer
    to the scaler to allow custom behaviour for single/multiple datasets."""
    self.scaler.update_for_minimisation(self.apm, self.curvatures)

  def get_num_matches(self):
    """Return the number of reflections in current minimisation cycle."""
    return self.scaler.Ih_table.size

  def rmsds(self):
    """Calculate RMSDs for the matches."""
    if self.scaler.Ih_table.blocked_Ih_table:
      R = flex.double([])
      n = 0
      for block in self.scaler.Ih_table.blocked_Ih_table.block_list:
        R.extend((self.calculate_residuals(block)**2) * block.weights)
        n += block.size
      restr = self.scaler.calculate_restraints(self.apm)
      #restr = self.scaler.compute_restraints_residuals_jacobian(self.apm)
      print(flex.sum(R))
      if restr:
        R.extend(restr[0])
      print(flex.sum(R))
      print(n)
      self._rmsds = [(flex.sum((R))/n)**0.5]
    else:
      R = (self.calculate_residuals(self.scaler.Ih_table)**2) * self.weights
      restr = self.scaler.calculate_restraints(self.apm)
      print(flex.sum(R))
      if restr:
        R.extend(restr[0])
      print(flex.sum(R))
      print(self.scaler.Ih_table.size)
      self._rmsds = [(flex.sum((R))/self.scaler.Ih_table.size)**0.5]
    if self.scaler.params.scaling_options.use_free_set:
      self.calculate_free_rmsds()
    return self._rmsds

  @staticmethod
  def achieved():
    """Method required by refinement engine."""
    return False #implement a method here?

  def calculate_residuals(self, Ih_table):
    """Return the residual vector."""
    R = Ih_table.intensities - (Ih_table.inverse_scale_factors
      * Ih_table.Ih_values)
    return R

  def calculate_free_rmsds(self):
    """Calculate an RMSD from the free set."""
    Ih_tab = self.scaler.Ih_table.free_Ih_table
    R = Ih_tab.intensities - (Ih_tab.inverse_scale_factors * Ih_tab.Ih_values)
    self._rmsds.append((flex.sum((R**2) * Ih_tab.weights)/Ih_tab.size)**0.5)

  def calculate_gradients(self, Ih_table):
    """Return a gradient vector on length len(self.apm.x)."""
    gsq = ((Ih_table.inverse_scale_factors)**2) * Ih_table.weights
    sumgsq = gsq * Ih_table.h_index_matrix
    rhl = Ih_table.intensities - (Ih_table.Ih_values
      * Ih_table.inverse_scale_factors)
    dIh = ((Ih_table.intensities - (Ih_table.Ih_values * 2.0 *
      Ih_table.inverse_scale_factors)) * Ih_table.weights)
    dIh_g = row_multiply(Ih_table.derivatives, dIh)
    dIh_red = dIh_g.transpose() * Ih_table.h_index_matrix
    dIh_by_dpi = row_multiply(dIh_red.transpose(), 1.0/sumgsq)
    term_1 = (-2.0 * rhl * Ih_table.weights * Ih_table.Ih_values *
              Ih_table.derivatives)
    term_2 = (-2.0 * rhl * Ih_table.weights * Ih_table.inverse_scale_factors *
              Ih_table.h_index_matrix) * dIh_by_dpi
    gradient = term_1 + term_2
    return gradient

  def calculate_jacobian(self, Ih_table):
    """Calculate the jacobian matrix, size Ih_table.size by len(self.apm.x)."""
    invsumgsq = 1.0 / ((Ih_table.weights * ((Ih_table.inverse_scale_factors)**2))
      * Ih_table.h_index_matrix)
    dIh = ((Ih_table.intensities - (Ih_table.Ih_values * 2.0 *
      Ih_table.inverse_scale_factors))) * Ih_table.weights
    dIh_g = row_multiply(Ih_table.derivatives, dIh)
    dIh_red = dIh_g.transpose() * Ih_table.h_index_matrix
    dIh_by_dpi = row_multiply(dIh_red.transpose(), invsumgsq)
    dIh_by_dpi = dIh_by_dpi.transpose() * Ih_table.h_expand_matrix
    term1 = row_multiply(dIh_by_dpi.transpose(),
      -1.0 * Ih_table.inverse_scale_factors)
    term2 = row_multiply(Ih_table.derivatives, Ih_table.Ih_values)
    for i, col in enumerate(term2.cols()):
      term1[:, i] -= col # Sum columns to create the jacobian
    return term1 # Return the jacobian

  def calculate_curvatures(self, Ih_table):
    """Return the second derivative of the target function."""
    assert 0, "method not yet implemented."
    """Calculate the second derivatives"""
    # Make some shorthand notation
    '''Ih_tab = self.scaler.Ih_table
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
    Ihprime = row_multiply(Ihprime_numerator, v_inv)
    #n_unique_groups x n_params matrix

    Ihprime_expanded = (Ihprime.transpose() * Ih_tab.h_expand_matrix).transpose()
    a = row_multiply(gprime, Ih)
    b = row_multiply(Ihprime_expanded, g)
    for j, col in enumerate(a.cols()):
      b[:, j] += col #n_refl x n_params
    rprime = b
    rprimesq = elementwise_square(rprime)
    # First term is sum_refl 2.0 * (r')2
    first_term = 2.0 * w * rprimesq
    #print(list(first_term))

    # Second term is sum_refl -2.0 * r * r''
    #  = -2.0 * r * (Ih g'' + 2g'Ih' + g Ih'') = A + B + C
    # Ih '' = u'/v - Ih v'/v. so C = -2 r g u'/v + 2 r g Ih' v'/v = D + E
    r = I - (g * Ih)
    if self.apm.curvatures:
      A = -2.0 * w * r * Ih * self.apm.curvatures #len n_params
    B = row_multiply(gprime, (-4.0 * w * r)).transpose() * Ih_tab.h_index_matrix
    B = B.transpose() * Ihprime_expanded #problem -
    # need to do elementwise multiplication of two big matrices.

    reduced_prefactor = (-2.0 * r * g * Ih_tab.h_index_matrix)
    B = -4.0 * r * ((h_idx_transpose * gprime).transpose() * Ihprime) #len n_params
    vprime_over_v = row_multiply((h_idx_transpose * row_multiply(gprime,
      (2.0 * w * g))), v_inv) #len n_unique_groups x n_params
    E = (-1.0 * reduced_prefactor * Ihprime) * vprime_over_v

    #E = (2.0 * r * g * Ih_tab.h_index_matrix) * row_multiply(Ihprime,
    # vprime_over_v) #len n_params

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
      return curvs
    curvs = B + D2 + D3 + E + first_term

    return curvs'''

  # The following methods are for adaptlbfgs.
  def compute_functional_gradients(self, block=None):
    """Return the functional and gradients."""
    if block:
      # block is an Ih_table block
      resids = self.calculate_residuals(block)
      gradients = self.calculate_gradients(block)
      weights = block.weights
    else:
      resids = self.calculate_residuals(self.scaler.Ih_table)
      gradients = self.calculate_gradients(self.scaler.Ih_table)
      weights = self.scaler.Ih_table.weights
    functional = flex.sum(resids**2 * weights)
    return functional, gradients

  def compute_functional_gradients_and_curvatures(self, block=None):
    """Return the functional, gradients and curvatures."""
    if block:
      resids = self.calculate_residuals(block)
      gradients = self.calculate_gradients(block)
      curvatures = self.calculate_curvatures(block)
      weights = block.weights
    else:
      resids = self.calculate_residuals(self.scaler.Ih_table)
      gradients = self.calculate_gradients(self.scaler.Ih_table)
      curvatures = self.calculate_curvatures(self.scaler.Ih_table)
      weights = self.scaler.Ih_table.weights
    functional = flex.sum(resids**2 * weights)
    return functional, gradients, curvatures

  def compute_restraints_functional_gradients_and_curvatures(self):
    """Return the restrains for functional, gradients and curvatures."""
    restraints = None
    #if 'absorption' in self.apm.components_list:
    restr = self.scaler.calculate_restraints(self.apm)
    if restr:
      resid_restr = flex.sum(restr[0]) #want just a value to add to total functional here
      grad_restr = restr[1]
      restraints = [resid_restr, grad_restr, None]
    return restraints #list of restraints to add to resid, grads and curvs?

  # The following methods are for adaptlstbx (GN/ LM algorithms)
  def compute_residuals(self, block=None):
    """Return the residuals array and weights."""
    if block:
      residuals = self.calculate_residuals(block)
      weights = block.weights
    else:
      residuals = self.calculate_residuals(self.scaler.Ih_table)
      weights = self.weights
    return residuals, weights

  def compute_residuals_and_gradients(self, block=None):
    """Return the residuals array, jacobian matrix and weights."""
    if block:
      residuals = self.calculate_residuals(block)
      jacobian = self.calculate_jacobian(block)
      weights = block.weights
    else:
      residuals = self.calculate_residuals(self.scaler.Ih_table)
      jacobian = self.calculate_jacobian(self.scaler.Ih_table)
      weights = self.weights
    return residuals, jacobian, weights

  def compute_restraints_residuals_and_gradients(self):
    """Return the restraints for the residuals and jacobian."""
    restraints = None
    if 'absorption' in self.apm.components_list:
      restraints = self.scaler.compute_restraints_residuals_jacobian(self.apm)
    return restraints

class ScalingTargetFixedIH(ScalingTarget):
  """A special implementation of scaling target for when the scaling is to be
  done against a fixed reference Ih set (i.e scaler is a TargetScaler)
  """
  _grad_names = ['dI_dp']
  rmsd_names = ["RMSD_I"]
  rmsd_units = ["a.u"]

  def __init__(self, scaler, apm, curvatures=False):
    self.scaler = scaler
    self.apm = apm
    self.weights = flex.double([])
    for scaler in self.scaler.unscaled_scalers:
      self.weights.extend(scaler.Ih_table.weights)
    self.curvatures = curvatures
    if self.scaler.params.scaling_options.use_free_set:
      self.rmsd_names = ["RMSD_I", "Free RMSD_I"]
      self.rmsd_units = ["a.u", "a.u"]
    # Quantities to cache each step
    self._rmsds = None

  def calculate_free_rmsds(self):
    """Calculate an RMSD from the free set."""
    RMSD = flex.double([])
    for unscaled_scaler in self.scaler.unscaled_scalers:
      Ih_tab = unscaled_scaler.Ih_table.free_Ih_table
      R = Ih_tab.intensities - (Ih_tab.inverse_scale_factors * Ih_tab.Ih_values)
      RMSD.extend((R**2) * Ih_tab.weights)
    self._rmsds.append((flex.sum(RMSD)/RMSD.size())**0.5)

  '''def calculate_gradients(self):
    G = flex.double([])
    for i, unscaled_scaler in enumerate(self.scaler.unscaled_scalers):
      apm = self.apm.apm_list[i]
      Ih_tab = unscaled_scaler.Ih_table
      rhl = Ih_tab.intensities - (Ih_tab.Ih_values * Ih_tab.inverse_scale_factors)
      G.extend(-2.0 * rhl * Ih_tab.weights * Ih_tab.Ih_values * apm.derivatives)
    return G'''

  def calculate_gradients(self, Ih_table):
    rhl = Ih_table.intensities - (Ih_table.Ih_values * Ih_table.inverse_scale_factors)
    G = -2.0 * rhl * Ih_table.weights * Ih_table.Ih_values * Ih_table.derivatives
    return G

  '''def calculate_residuals(self):
    R = flex.double([])
    for unscaled_scaler in self.scaler.unscaled_scalers:
      Ih_tab = unscaled_scaler.Ih_table
      R.extend(Ih_tab.intensities - (Ih_tab.inverse_scale_factors *
        Ih_tab.Ih_values))
    return R'''

  def calculate_jacobian(self):
    """Calculate the jacobian matrix, size Ih_table.size by len(self.apm.x)."""
    n_refl = sum([i.Ih_table.size for i in self.scaler.unscaled_scalers])
    jacobian = sparse.matrix(n_refl, self.apm.n_active_params)
    n_cumul_rows = 0
    for i, scaler in enumerate(self.scaler.unscaled_scalers):
      apm = self.apm.apm_list[i]
      j = row_multiply(apm.derivatives, -1.0 * scaler.Ih_table.Ih_values)
      jacobian.assign_block(j, n_cumul_rows, self.apm.apm_data[i]['start_idx'])
      n_cumul_rows += apm.n_obs
    return jacobian

  def compute_functional_gradients(self, block=None):
    """Return the functional and gradients."""
    if block:
      # block is an Ih_table block
      resids = self.calculate_residuals(block)
      gradients = self.calculate_gradients(block)
      weights = block.weights
    else:
      residuals = flex.double()
      weights = flex.double()
      gradients = flex.double()
      for unscaled_scaler in self.scaler.unscaled_scalers:
        residuals.extend(self.calculate_residuals(unscaled_scaler.Ih_table))
        gradients.extend(self.calculate_gradients(unscaled_scaler.Ih_table))
        weights.extend(unscaled_scaler.Ih_table.weights)
    functional = flex.sum(residuals**2 * weights)
    return functional, gradients

  def compute_functional_gradients_and_curvatures(self, block=None):
    """Return the functional, gradients and curvatures."""
    if block:
      resids = self.calculate_residuals(block)
      gradients = self.calculate_gradients(block)
      curvatures = self.calculate_curvatures(block)
      weights = block.weights
    else:
      residuals = flex.double()
      weights = flex.double()
      gradients = flex.double()
      curvatures = flex.double()
      for unscaled_scaler in self.scaler.unscaled_scalers:
        residuals.extend(self.calculate_residuals(unscaled_scaler.Ih_table))
        gradients.extend(self.calculate_gradients(unscaled_scaler.Ih_table))
        curvatures.extend(self.calculate_curvatures(unscaled_scaler.Ih_table))
        weights.extend(unscaled_scaler.Ih_table.weights)
    functional = flex.sum(residuals**2 * weights)
    return functional, gradients, curvatures

  def compute_restraints_functional_gradients_and_curvatures(self):
    """Return the restrains for functional, gradients and curvatures."""
    restraints = None
    #if 'absorption' in self.apm.components_list:
    restr = self.scaler.calculate_restraints(self.apm)
    if restr:
      resid_restr = flex.sum(restr[0]) #want just a value to add to total functional here
      grad_restr = restr[1]
      restraints = [resid_restr, grad_restr, None]
    return restraints #list of restraints to add to resid, grads and curvs?

  # The following methods are for adaptlstbx (GN/ LM algorithms)
  def compute_residuals(self, block=None):
    """Return the residuals array and weights."""
    if block:
      residuals = self.calculate_residuals(block)
      weights = block.weights
    else:
      residuals = flex.double()
      weights = flex.double()
      for unscaled_scaler in self.scaler.unscaled_scalers:
        residuals.extend(self.calculate_residuals(unscaled_scaler.Ih_table))
        weights.extend(unscaled_scaler.Ih_table.weights)
    return residuals, weights

  def compute_residuals_and_gradients(self, block=None):
    """Return the residuals array, jacobian matrix and weights."""
    if block:
      residuals = self.calculate_residuals(block)
      jacobian = self.calculate_jacobian(block)
      weights = block.weights
    else:
      residuals = self.calculate_residuals(self.scaler.Ih_table)
      jacobian = self.calculate_jacobian(self.scaler.Ih_table)
      weights = self.weights
    return residuals, jacobian, weights

  def compute_restraints_residuals_and_gradients(self):
    """Return the restraints for the residuals and jacobian."""
    restraints = None
    if 'absorption' in self.apm.components_list:
      restraints = self.scaler.compute_restraints_residuals_jacobian(self.apm)
    return restraints
