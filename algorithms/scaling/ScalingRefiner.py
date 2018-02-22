"""Contains classes for scaling refinement engines - inherited from the 
dials.refinement engine with a few methods overwritten to use them with scaling
code"""

from __future__ import absolute_import, division
import logging
from dials.algorithms.refinement.engine import SimpleLBFGS,\
 GaussNewtonIterations, LevenbergMarquardtIterations
#logger = logging.getLogger(__name__)
logger = logging.getLogger('dials')



TARGET_ACHIEVED = "RMSD target achieved"
RMSD_CONVERGED = "RMSD no longer decreasing"
STEP_TOO_SMALL = "Step too small"
OBJECTIVE_INCREASE = "Refinement failure: objective increased"
MAX_ITERATIONS = "Reached maximum number of iterations"
MAX_TRIAL_ITERATIONS = "Reached maximum number of consecutive unsuccessful trial steps"
DOF_TOO_LOW = "Not enough degrees of freedom to refine"

"""refinery_phil_str = '''
refinery
  .help = "Parameters to configure the refinery"
  .expert_level = 1
{
  engine = SimpleLBFGS LBFGScurvs GaussNewton *LevMar SparseLevMar
    .help = "The minimisation engine to use"
    .type = choice

  max_iterations = None
    .help = "Maximum number of iterations in refinement before termination."
            "None implies the engine supplies its own default."
    .type = int(value_min=1)

  log = None
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
scaling_refinery_phil_scope = parse(refinery_phil_str)"""


class ScalingRefinery(object):
  'mixin class to add extra return method'
  def __init__(self, scaler):
    self._scaler = scaler

  def print_step_table(self):
    """print useful output about refinement steps in the form of a simple table"""

    from libtbx.table_utils import simple_table

    logger.info("\nRefinement steps:")

    header = ["Step", "Nref"]
    for (name, units) in zip(self._target.rmsd_names, self._target.rmsd_units):
      header.append(name + "\n(" + units + ")")

    rows = []
    for i in range(self.history.get_nrows()):
      rmsds = [r for r in self.history["rmsd"][i]]
      rows.append([str(i), str(self.history["num_reflections"][i])] + \
        ["%.5g" % r for r in rmsds])

    st = simple_table(rows, header)
    logger.info(st.format())
    logger.info(self.history.reason_for_termination)


  def return_scaler(self):
    '''return scaler method'''
    from dials.algorithms.scaling.Scaler import MultiScalerBase, SingleScalerBase
    self.print_step_table()

    if isinstance(self._scaler, SingleScalerBase):
      if self._parameters.var_cov_matrix:
        self._scaler.update_var_cov(self._parameters)
        #self._scaler.var_cov_matrix = self._parameters.var_cov_matrix
    elif self._scaler.id_ == 'multi':
      if self._parameters.apm_list[0].var_cov_matrix:
        for i, scaler in enumerate(self._scaler.single_scalers):
          scaler.update_var_cov(self._parameters.apm_list[i])
          #scaler.var_cov_matrix = self._parameters.apm_list[i].var_cov_matrix
    elif self._scaler.id_ == 'target':
      if self._parameters.apm_list[0].var_cov_matrix:
        for i, scaler in enumerate(self._scaler.unscaled_scalers):
          scaler.update_var_cov(self._parameters.apm_list[i])
          #scaler.var_cov_matrix = self._parameters.apm_list[i].var_cov_matrix

    if not isinstance(self._scaler, MultiScalerBase):
      if 'scale' in self._parameters.components:
        self._scaler.normalise_scale_component()
      if 'decay' in self._parameters.components:
        self._scaler.normalise_decay_component()
    return self._scaler


class ScalingSimpleLBFGS(SimpleLBFGS, ScalingRefinery):
  """Adapt Refinery for L-BFGS minimiser"""
  def __init__(self, scaler, *args, **kwargs):
    super(ScalingSimpleLBFGS, self).__init__(*args, **kwargs)
    ScalingRefinery.__init__(self, scaler)

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
        #print(flex.sum(residuals)+flex.sum(restraints[0]))
    return

class ScalingGaussNewtonIterations(ScalingLstbxBuildUpMixin, GaussNewtonIterations):
  """Refinery implementation, using lstbx Gauss Newton iterations"""

  # defaults that may be overridden
  gradient_threshold = 1.e-10
  step_threshold = None
  damping_value = 0.0007
  max_shift_over_esd = 15
  convergence_as_shift_over_esd = 1e-5

  def __init__(self, scaler, target, prediction_parameterisation, constraints_manager=None,
               log=None, verbosity=0, tracking=None,
               max_iterations=20, **kwds):

    GaussNewtonIterations.__init__(
             self, target, prediction_parameterisation, constraints_manager,
             log=log, verbosity=verbosity, tracking=tracking,
             max_iterations=max_iterations)
    ScalingLstbxBuildUpMixin.__init__(self, scaler)

class ScalingLevenbergMarquardtIterations(ScalingLstbxBuildUpMixin, LevenbergMarquardtIterations):
  """Refinery implementation, employing lstbx Levenberg Marquadt
  iterations"""

  def __init__(self, scaler, target, prediction_parameterisation, constraints_manager=None,
               log=None, verbosity=0, tracking=None,
               max_iterations=20, **kwds):

    LevenbergMarquardtIterations.__init__(
             self, target, prediction_parameterisation, constraints_manager,
             log=log, verbosity=verbosity, tracking=tracking,
             max_iterations=max_iterations)
    ScalingLstbxBuildUpMixin.__init__(self, scaler)
