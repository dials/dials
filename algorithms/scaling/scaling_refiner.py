""" Classes for scaling refinement engines.

Classes are inherited from the dials.refinement engine with a few
methods overwritten to use them with scaling code."""

from __future__ import absolute_import, division
import logging
from dials.algorithms.refinement.engine import SimpleLBFGS,\
 GaussNewtonIterations, LevenbergMarquardtIterations, LBFGScurvs
from libtbx.phil import parse
#logger = logging.getLogger(__name__)
logger = logging.getLogger('dials')


TARGET_ACHIEVED = "RMSD target achieved"
RMSD_CONVERGED = "RMSD no longer decreasing"
STEP_TOO_SMALL = "Step too small"
OBJECTIVE_INCREASE = "Refinement failure: objective increased"
MAX_ITERATIONS = "Reached maximum number of iterations"
MAX_TRIAL_ITERATIONS = "Reached maximum number of consecutive unsuccessful trial steps"
DOF_TOO_LOW = "Not enough degrees of freedom to refine"

scaling_refinery_phil_str = '''
scaling_refinery
  .help = "Parameters to configure the refinery"
  .expert_level = 1
{
  engine = *SimpleLBFGS GaussNewton LevMar
    .help = "The minimisation engine to use for the main scaling algorithm"
    .type = choice
  max_iterations = None
    .help = "Maximum number of iterations in refinement before termination."
            "None implies the engine supplies its own default."
    .type = int(value_min=1)

  full_matrix_engine = GaussNewton *LevMar
    .help = "The minimisation engine to use for a full matrix round of
             minimisation after the main scaling, in order to determine
             error estimates."
    .type = choice
  full_matrix_max_iterations = None
    .help = "Maximum number of iterations before termination in the full matrix
             minimisation round."
            "None implies the engine supplies its own default."
    .type = int(value_min=1)
}
'''
'''refinery_log = None
    .help = "Filename for an optional log that a minimisation engine may use"
            "to write additional information"
    .type = path

  journal
    .help = "Extra items to track in the refinement history"
  {
    track_step = False
      .help = "Record parameter shifts history in the refinement journal, if"
              "the engine supports it."
      .type = bool

    track_gradient = False
      .help = "Record parameter gradients history in the refinement journal, if"
              "the engine supports it."
      .type = bool

    track_parameter_correlation = False
      .help = "Record correlation matrix between columns of the Jacobian for"
              "each step of refinement."
      .type = bool

    track_condition_number = False
      .help = "Record condition number of the Jacobian for each step of "
              "refinement."
      .type = bool

    track_out_of_sample_rmsd = False
      .type = bool
      .help = "Record RMSDs calculated using the refined experiments with"
              "reflections not used in refinement at each step. Only valid if a"
              "subset of input reflections was taken for refinement"
  }
}
'''
scaling_refinery_phil_scope = parse(scaling_refinery_phil_str)

def scaling_refinery(engine, target, prediction_parameterisation,
    max_iterations):
  """Return the correct engine based on phil parameters."""
  if engine == 'SimpleLBFGS':
    return ScalingSimpleLBFGS(target=target,
      prediction_parameterisation=prediction_parameterisation,
      max_iterations=max_iterations)
  elif engine == 'LBFGScurvs':
    return ScalingLBFGScurvs(target=target,
      prediction_parameterisation=prediction_parameterisation,
      max_iterations=max_iterations)
  elif engine == 'GaussNewton':
    return ScalingGaussNewtonIterations(target=target,
      prediction_parameterisation=prediction_parameterisation,
      max_iterations=max_iterations)
  elif engine == 'LevMar':
    return ScalingLevenbergMarquardtIterations(target=target,
      prediction_parameterisation=prediction_parameterisation,
      max_iterations=max_iterations)

def error_model_refinery(engine, target, max_iterations):
  """Return the correct engine based on phil parameters.

  Note that here the target also takes the role of the predication
  parameterisation by implementing the set_param_vals and get_param_vals
  methods (the code is organised in this way to allow the use of the
  dials.refinement engines)."""
  if engine == 'SimpleLBFGS':
    return ErrorModelSimpleLBFGS(target=target,
      prediction_parameterisation=target,
      max_iterations=max_iterations)
  '''elif engine == 'LBFGScurvs':
    assert 0, 'LBFGS with curvatures not yet implemented'
  elif engine == 'GaussNewton':
    return ScalingGaussNewtonIterations(target=target,
      prediction_parameterisation=prediction_parameterisation,
      max_iterations=max_iterations)
  elif engine == 'LevMar':
    return ScalingLevenbergMarquardtIterations(target=target,
      prediction_parameterisation=prediction_parameterisation,
      max_iterations=max_iterations)'''

def print_step_table(refinery):
  """print useful output about refinement steps in the form of a simple table"""

  from libtbx.table_utils import simple_table

  logger.info("\nRefinement steps:")

  header = ["Step", "Nref"]
  for (name, units) in zip(refinery._target.rmsd_names, refinery._target.rmsd_units):
    header.append(name + "\n(" + units + ")")

  rows = []
  for i in range(refinery.history.get_nrows()):
    rmsds = [r for r in refinery.history["rmsd"][i]]
    rows.append([str(i), str(refinery.history["num_reflections"][i])] + \
      ["%.5g" % r for r in rmsds])

  st = simple_table(rows, header)
  logger.info(st.format())
  logger.info(refinery.history.reason_for_termination)

class ErrorModelRefinery(object):
  """Mixin class to add extra return method."""
  def __init__(self, error_model_manager):
    self.error_manager = error_model_manager

  def return_error_manager(self):
    """Set error manager parameters and return error manager."""
    print_step_table(self)
    self.error_manager.refined_parameters = self._target.x
    logger.info("\nMinimised error model with parameters {0:.5f} and {1:.5f}. {sep}"
          .format(self._target.x[0], abs(self._target.x[1]), sep='\n'))
    return self.error_manager


class ScalingRefinery(object):
  'mixin class to add extra return method'
  def __init__(self, scaler):
    self._scaler = scaler

  def return_scaler(self):
    '''return scaler method'''
    from dials.algorithms.scaling.scaler import MultiScalerBase
    print_step_table(self)

    if self._scaler.id_ == 'single':
      if self._parameters.var_cov_matrix:
        self._scaler.update_var_cov(self._parameters)
        self._scaler.experiments.scaling_model.set_scaling_model_as_scaled()
    elif self._scaler.id_ == 'multi' or self._scaler.id_ == 'target':
      if self._parameters.apm_list[0].var_cov_matrix: #test if has been set
        for i, scaler in enumerate(self._scaler.single_scalers):
          scaler.update_var_cov(self._parameters.apm_list[i])
          scaler.experiments.scaling_model.set_scaling_model_as_scaled()

    if not isinstance(self._scaler, MultiScalerBase):
      self._scaler.experiments.scaling_model.normalise_components()
    return self._scaler


class ScalingSimpleLBFGS(SimpleLBFGS, ScalingRefinery):
  """Adapt Refinery for L-BFGS minimiser"""
  def __init__(self, *args, **kwargs):
    super(ScalingSimpleLBFGS, self).__init__(*args, **kwargs)
    ScalingRefinery.__init__(self, self._target.scaler)

  def compute_functional_gradients_and_curvatures(self):
    """overwrite method to avoid calls to 'blocks' methods of target"""
    self.prepare_for_step()

    f, g = self._target.compute_functional_gradients()

    # restraints terms
    restraints = \
      self._target.compute_restraints_functional_gradients_and_curvatures()

    if restraints:
      f += restraints[0]
      g += restraints[1]

    return f, g, None

class ScalingLBFGScurvs(LBFGScurvs, ScalingRefinery):
  """Adapt Refinery for L-BFGS minimiser"""
  def __init__(self, *args, **kwargs):
    super(ScalingLBFGScurvs, self).__init__(*args, **kwargs)
    ScalingRefinery.__init__(self, self._target.scaler)
    self._target.curvatures = True

  def compute_functional_gradients_and_curvatures(self):
    """overwrite method to avoid calls to 'blocks' methods of target"""
    self.prepare_for_step()

    f, g, curv = self._target.compute_functional_gradients_and_curvatures()

    # restraints terms
    restraints = \
      self._target.compute_restraints_functional_gradients_and_curvatures()

    if restraints:
      f += restraints[0]
      g += restraints[1]

    return f, g, curv


class ErrorModelSimpleLBFGS(SimpleLBFGS, ErrorModelRefinery):
  """Adapt Refinery for L-BFGS minimiser"""
  def __init__(self, *args, **kwargs):
    super(ErrorModelSimpleLBFGS, self).__init__(*args, **kwargs)
    ErrorModelRefinery.__init__(self, self)

  def compute_functional_gradients_and_curvatures(self):
    """overwrite method to avoid calls to 'blocks' methods of target"""
    self.prepare_for_step()

    f, g, _ = self._target.compute_functional_gradients_and_curvatures()

    # restraints terms
    restraints = \
      self._target.compute_restraints_functional_gradients_and_curvatures()

    if restraints:
      f += restraints[0]
      g += restraints[1]

    return f, g, None

class ScalingLstbxBuildUpMixin(ScalingRefinery):
  '''Mixin class to overwrite the build_up method in AdaptLstbx'''
  def build_up(self, objective_only=False):
    'overwrite method from Adaptlstbx'
    # set current parameter values
    self.prepare_for_step()

    # Reset the state to construction time, i.e. no equations accumulated
    self.reset()

    # observation terms
    if objective_only:
      residuals, weights = self._target.compute_residuals()
      self.add_residuals(residuals, weights)
    else:
      residuals, self._jacobian, weights = \
            self._target.compute_residuals_and_gradients()
      self.add_equations(residuals, self._jacobian, weights)

    restraints = self._target.compute_restraints_residuals_and_gradients()
    if restraints:
      if objective_only:
        self.add_residuals(restraints[0], restraints[2])
      else:
        self.add_equations(restraints[0], restraints[1], restraints[2])

    return

class ScalingGaussNewtonIterations(ScalingLstbxBuildUpMixin, GaussNewtonIterations):
  """Refinery implementation, using lstbx Gauss Newton iterations"""

  # defaults that may be overridden
  gradient_threshold = 1.e-10
  step_threshold = None
  damping_value = 0.0007
  max_shift_over_esd = 15
  convergence_as_shift_over_esd = 1e-5

  def __init__(self, target, prediction_parameterisation, constraints_manager=None,
               log=None, verbosity=0, tracking=None,
               max_iterations=20):

    GaussNewtonIterations.__init__(
             self, target, prediction_parameterisation, constraints_manager,
             log=log, verbosity=verbosity, tracking=tracking,
             max_iterations=max_iterations)
    ScalingLstbxBuildUpMixin.__init__(self, self._target.scaler)

class ScalingLevenbergMarquardtIterations(ScalingLstbxBuildUpMixin, LevenbergMarquardtIterations):
  """Refinery implementation, employing lstbx Levenberg Marquadt
  iterations"""

  def __init__(self, target, prediction_parameterisation, constraints_manager=None,
               log=None, verbosity=0, tracking=None,
               max_iterations=20):

    LevenbergMarquardtIterations.__init__(
             self, target, prediction_parameterisation, constraints_manager,
             log=log, verbosity=verbosity, tracking=tracking,
             max_iterations=max_iterations)
    ScalingLstbxBuildUpMixin.__init__(self, self._target.scaler)
