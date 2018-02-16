"""Contains classes for refinement engines. Refinery is the shared interface,
LevenbergMarquardtIterations, GaussNewtonIterations, SimpleLBFGS and LBFGScurvs
are the current concrete implementations"""

from __future__ import absolute_import, division
import logging
#logger = logging.getLogger(__name__)
logger = logging.getLogger('dials')

from scitbx import lbfgs
from scitbx.array_family import flex
import libtbx
from libtbx import easy_mp
from libtbx.phil import parse

# use lstbx classes
from scitbx.lstbx import normal_eqns, normal_eqns_solving

class Refinery(object):

  def __init__(self, scaler, apm, target=None,
      constraints_manager=None, max_iterations = None):
    self._scaler = scaler
    self._apm = apm
    self._target = target
    self._constr_manager = constraints_manager
    
    self.x = self._apm.x
    if self._constr_manager is not None:
      self.x = self._constr_manager.constrain_parameters(self.x)

    self._f = None
    self._g = None
    self._jacobian = None

    self._target_achieved = False
    self._max_iterations = max_iterations

    self.prepare_for_step()

  def prepare_for_step(self):
    '( i.e. update_for_minimisation)'
    """Update the parameterisation and prepare the target function"""
    self._scaler.update_for_minimisation(self._apm)

  def test_for_termination(self):
    """Return True if refinement should be terminated"""

    # Basic version delegate to the Target class. Derived classes may
    # implement other termination criteria
    self._target_achieved = self._target.achieved()

    return self._target_achieved

  def run(self):
    """
    To be implemented by derived class. It is expected that each step of
    refinement be preceeded by a call to prepare_for_step and followed by
    calls to update_journal and test_for_termination (in that order).
    """
    # Specify a minimizer and its parameters, and run
    raise NotImplementedError()

  def return_scaler(self):
    '''return scaler method'''
    from dials.algorithms.scaling.Scaler import MultiScalerBase
    if not isinstance(self._scaler, MultiScalerBase):
      if 'scale' in self._apm.components:
        self._scaler.normalise_scale_component()
      if 'decay' in self._apm.components:
        self._scaler.normalise_decay_component()
    return self._scaler


class AdaptLbfgs(Refinery):
  """Adapt Refinery for L-BFGS minimiser"""
  def __init__(self, *args, **kwargs):
    super(AdaptLbfgs, self).__init__(*args, **kwargs)

    self._termination_params = lbfgs.termination_parameters(
      max_iterations = self._max_iterations)

  def compute_functional_and_gradients(self):
    self.prepare_for_step()

    L, dL_dp, = self._scaler.get_target_function(self._apm)
    self._f = L
    self._g = dL_dp
    logger.info("Residual sum: %12.6g" % self._f)
    return self._f, self._g


  def run(self, curvatures=False):
    """
    Run the minimiser, keeping track of its log.
    """
    if curvatures: self.diag_mode = "always"
    self.minimizer = lbfgs.run(target_evaluator=self,
        termination_params=self._termination_params)

class AdaptLstbx(Refinery, normal_eqns.non_linear_ls,
    normal_eqns.non_linear_ls_mixin):
  """Adapt Refinery for lstbx"""

  def __init__(self, scaler, apm, target=None,
      constraints_manager=None, max_iterations=None):
    super(AdaptLstbx, self).__init__(scaler, apm, target=target,
      constraints_manager=constraints_manager, max_iterations=max_iterations)

    # required for restart to work (do I need that method?)
    self.x_0 = self.x.deep_copy()

    # keep attribute for the Cholesky factor required for ESD calculation
    self.cf = None

    normal_eqns.non_linear_ls.__init__(self, n_parameters = len(self.x))

  def restart(self):
    self.x = self.x_0.deep_copy()
    self.old_x = None

  def parameter_vector_norm(self):
    return self.x.norm()

  def build_up(self, objective_only=False):

    # code here to calculate the residuals. Rely on the target class
    # for this

    # I need to use the weights. They are the variances of the
    # observations... See http://en.wikipedia.org/wiki/Non-linear_least_squares
    # at 'diagonal weight matrix'

    # set current parameter values
    self.prepare_for_step()

    # Reset the state to construction time, i.e. no equations accumulated
    self.reset()

    # observation terms
    if objective_only:
      residual, _, = self._scaler.get_target_function(self._apm)
      residuals, weights = self._target.compute_residuals()
      self.add_residuals(residuals, weights)
    else:
      residuals, self._jacobian, weights = self._scaler.get_residuals_jacobian_weight(self._apm)
      print(flex.sum(residuals))
        #self._target.compute_residuals_and_gradients(block)
      j = self._jacobian
      #if self._constr_manager is not None:
      #  j = self._constr_manager.constrain_jacobian(j)
      self.add_equations(residuals, j, weights)

  def step_forward(self):
    self.old_x = self.x.deep_copy()
    self.x += self.step()

  def step_backward(self):
    if self.old_x is None:
      return False
    else:
      self.x, self.old_x = self.old_x, None
      return True

  def set_cholesky_factor(self):
    """Set the Cholesky factor required for ESD calculation. This method is
    valid only for the LSTBX dense matrix interface"""

    self.cf = self.step_equations().cholesky_factor_packed_u().deep_copy()

  '''def calculate_esds(self):
    """Calculate ESDs of parameters"""

    # it is possible to get here with zero steps taken by the minimiser. For
    # example by failing for the MAX_TRIAL_ITERATIONS reason before any forward
    # steps are taken with the LevMar engine. If so the below is invalid,
    # so return early
    if self.history.get_nrows() == 0: return None

    if self.cf is None: return None

    # if constraints were used then the normal matrix has fewer rows/columns
    # than the number of expanded parameters. At the moment, do not support
    # this calculation when constraints were used
    if self._constr_manager is not None: return None

    # invert normal matrix from N^-1 = (U^-1)(U^-1)^T
    cf_inv = self.cf.matrix_packed_u_as_upper_triangle().\
        matrix_inversion()
    nm_inv = cf_inv.matrix_multiply_transpose(cf_inv)

    # keep the estimated parameter variance-covariance matrix
    self.parameter_var_cov = \
        self.history["reduced_chi_squared"][-1] * nm_inv
    # send this back to the models to calculate their uncertainties
    self._parameters.calculate_model_state_uncertainties(
      self.parameter_var_cov)

    # send parameter variances back to the parameter classes
    # themselves, for reporting purposes and for building restraints
    # based on existing parameterisations.
    s2 = self.parameter_var_cov.matrix_diagonal()
    assert s2.all_ge(0.0)
    s = flex.sqrt(s2)
    self._parameters.set_param_esds(s)'''

  def _print_normal_matrix(self):
    """Print the full normal matrix at the current step. For debugging only"""
    logger.debug("The normal matrix for the current step is:")
    logger.debug(self.normal_matrix_packed_u().\
          matrix_packed_u_as_symmetric().\
          as_scitbx_matrix().matlab_form(format=None,
          one_row_per_line=True))
    logger.debug("\n")

class GaussNewtonIterations(AdaptLstbx, normal_eqns_solving.iterations):
  """Refinery implementation, using lstbx Gauss Newton iterations"""

  # defaults that may be overridden
  gradient_threshold = 1.e-10
  step_threshold = None
  damping_value = 0.0007
  max_shift_over_esd = 15
  convergence_as_shift_over_esd = 1e-5

  def __init__(self, scaler, apm, target=None,
      constraints_manager=None, max_iterations=20, **kwds):

    AdaptLstbx.__init__(self, scaler, apm, target=target,
      constraints_manager=constraints_manager,
      max_iterations=max_iterations)

    # add an attribute to the journal
    #self.history.add_column("reduced_chi_squared")#flex.double()

    # adopt any overrides of the defaults above
    #libtbx.adopt_optional_init_args(self, kwds)

  def run(self):
    self.n_iterations = 0

    # prepare for first step
    self.build_up()

    # return early if refinement is not possible
    if self.dof < 1:
      #self.history.reason_for_termination = DOF_TOO_LOW
      return

    while True:

      # set functional and gradients for the step (to add to the history)
      self._f = self.objective()
      self._g = -self.opposite_of_gradient()

      # cache some items for the journal prior to solve
      pvn = self.parameter_vector_norm()
      gn = self.opposite_of_gradient().norm_inf()

      # solve the normal equations
      self.solve()

      '''# standard journalling
      self.update_journal()
      logger.debug("Step %d", self.history.get_nrows() - 1)

      # add cached items to the journal
      self.history.set_last_cell("parameter_vector_norm", pvn)
      self.history.set_last_cell("gradient_norm", gn)

      # extra journalling post solve
      if "solution" in self.history:
        self.history.set_last_cell("solution", self.actual.step().deep_copy())
      self.history.set_last_cell("solution_norm", self.step().norm())
      self.history.set_last_cell("reduced_chi_squared", self.chi_sq())'''

      # test termination criteria
      '''if self.test_for_termination():
        #self.history.reason_for_termination = TARGET_ACHIEVED
        break

      if self.test_rmsd_convergence():
        #self.history.reason_for_termination = RMSD_CONVERGED
        break

      if self.had_too_small_a_step():
        #self.history.reason_for_termination = STEP_TOO_SMALL
        break'''

      '''if self.test_objective_increasing_but_not_nref():
        self.history.reason_for_termination = OBJECTIVE_INCREASE
        if self.step_backward():
          self.history.reason_for_termination += ". Parameters set back one step"
        self.prepare_for_step()
        break'''

      if self.n_iterations == self._max_iterations:
        #self.history.reason_for_termination = MAX_ITERATIONS
        break

      # prepare for next step
      self.step_forward()
      self.n_iterations += 1
      self.build_up()

    self.set_cholesky_factor()
    #self.calculate_esds()

    return

