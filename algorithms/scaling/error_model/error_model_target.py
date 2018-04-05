from dials.array_family import flex

class ErrorModelTarget(object):
  """Target function for a basic two parameter error model."""

  _grad_names = ['d(deltahl)_dp']
  rmsd_names = ["RMSD_deltahl"]
  rmsd_units = ["a.u"]

  def __init__(self, error_manager, starting_values=None):
    if not starting_values:
      starting_values = [1.0, 0.05]
    # Note - don't initialise with b = 0.0 or it gets stuck on 0.0!!
    self.error_model = error_manager
    self.x = starting_values

    # Quantities to cache each step
    self._rmsds = None

  def set_param_vals(self, x):
    '''method for refinement engine access'''
    self.x = x

  def get_param_vals(self):
    '''method for refinement engine access'''
    return self.x

  def predict(self):#do basis function calcuation and update Ih
    self.error_model.update_for_minimisation(self.x)

  def get_num_matches(self):
    return self.error_model.Ih_table.size

  def rmsds(self):
    """calculate unweighted RMSDs for the matches"""
    # cache rmsd calculation for achieved test
    R = self.calculate_residuals()
    R.extend(flex.double([self.compute_restraints_functional_gradients_and_curvatures()[0]]))
    self._rmsds = [(flex.sum((R))/self.error_model.Ih_table.size)**0.5]
    return self._rmsds

  def achieved(self):
    """Method required for refinement engine."""
    return False #implement a method here?

  def calculate_residuals(self):
    """Return the residual vector"""
    bin_vars = self.error_model.bin_variances
    R = ((flex.double(bin_vars.size(), 1.0) - bin_vars)**2)
    return R

  def calculate_gradients(self):
    'calculate the gradient vector'
    I_hl = self.error_model.Ih_table.intensities
    bin_vars = self.error_model.bin_variances
    sum_matrix = self.error_model.summation_matrix
    bin_counts = self.error_model.bin_counts
    dsig_da = self.error_model.sigmaprime/self.x[0]
    dsig_dc = self.x[1] * (I_hl**2) * (self.x[0]**2) / self.error_model.sigmaprime
    ddelta_dsigma = -1.0 * self.error_model.delta_hl / self.error_model.sigmaprime
    dsig_list = [ddelta_dsigma * dsig_da, ddelta_dsigma * dsig_dc]
    gradient = flex.double([])
    for deriv in dsig_list:
      term1 = 2.0 * self.error_model.delta_hl * deriv * sum_matrix
      term2a = self.error_model.delta_hl * sum_matrix
      term2b = deriv * sum_matrix
      grad = (-2.0 * (flex.double(bin_vars.size(), 1.0) - bin_vars)
        * ((term1 / bin_counts) - (2.0 * term2a * term2b / (bin_counts**2))))
      gradient.append(flex.sum(grad))
    return gradient

  # The following methods are for adaptlbfgs.
  def compute_functional_gradients(self):
    return flex.sum(self.calculate_residuals()), self.calculate_gradients()

  def compute_functional_gradients_and_curvatures(self):
    return flex.sum(self.calculate_residuals()), self.calculate_gradients(), None #curvatures

  def compute_restraints_functional_gradients_and_curvatures(self):
    R1 = 25.0
    R2 = 400.0
    residual_restraint = ((R1*((1.0 - self.x[0])**2))
      + (R2*((0.0001 - self.x[1])**2)))
    gradient_restraint = flex.double([-2.0 * R1 * (1.0 - self.x[0]),
      R2 * 2.0 * self.x[1]])
    return [residual_restraint, gradient_restraint, None]
