""" Classes for scaling refinement engines.

Classes are inherited from the dials.refinement engine with a few
methods overwritten to use them with scaling code."""


from __future__ import annotations

import logging

from libtbx.phil import parse

from dials.algorithms.refinement.engine import (
    GaussNewtonIterations,
    LevenbergMarquardtIterations,
    SimpleLBFGS,
)
from dials.algorithms.scaling.scaling_utilities import log_memory_usage
from dials.util import tabulate

logger = logging.getLogger("dials")


TARGET_ACHIEVED = "RMSD target achieved"
RMSD_CONVERGED = "RMSD no longer decreasing"
STEP_TOO_SMALL = "Step too small"
OBJECTIVE_INCREASE = "Refinement failure: objective increased"
MAX_ITERATIONS = "Reached maximum number of iterations"
MAX_TRIAL_ITERATIONS = "Reached maximum number of consecutive unsuccessful trial steps"
DOF_TOO_LOW = "Not enough degrees of freedom to refine"

scaling_refinery_phil_str = """
scaling_refinery
  .help = "Parameters to configure the refinery"
  .expert_level = 1
{
  engine = *SimpleLBFGS GaussNewton LevMar
    .help = "The minimisation engine to use for the main scaling algorithm"
    .type = choice
  refinement_order = *concurrent consecutive
    .type = choice
    .help = "Choice of whether to refine all model components concurrently, or"
             "in a consecutive order as allowed/defined by the scaling model."
    .expert_level = 2
  max_iterations = None
    .help = "Maximum number of iterations in refinement before termination."
            "None implies the engine supplies its own default."
    .type = int(value_min=1)
  rmsd_tolerance = 0.0001
    .type = float(value_min=1e-6)
    .help = "Tolerance at which to stop scaling refinement. This is a relative
            value, the convergence criterion is (rmsd[i] - rmsd[i-1])/rmsd[i] <
            rmsd_tolerance."
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
"""
"""refinery_log = None
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
"""
scaling_refinery_phil_scope = parse(scaling_refinery_phil_str)


def scaling_refinery(
    engine, scaler, target, prediction_parameterisation, max_iterations
):
    """Return the correct engine based on phil parameters."""
    if engine == "SimpleLBFGS":
        return ScalingSimpleLBFGS(
            scaler,
            target=target,
            prediction_parameterisation=prediction_parameterisation,
            max_iterations=max_iterations,
        )
    elif engine == "GaussNewton":
        return ScalingGaussNewtonIterations(
            scaler,
            target=target,
            prediction_parameterisation=prediction_parameterisation,
            max_iterations=max_iterations,
        )
    elif engine == "LevMar":
        return ScalingLevenbergMarquardtIterations(
            scaler,
            target=target,
            prediction_parameterisation=prediction_parameterisation,
            max_iterations=max_iterations,
        )


def print_step_table(refinery):
    """print useful output about refinement steps in the form of a simple table"""

    logger.info("\nRefinement steps:")

    header = ["Step", "Nref"]
    for (name, units) in zip(refinery._target.rmsd_names, refinery._target.rmsd_units):
        header.append(name + "\n(" + units + ")")

    rows = []
    for i in range(refinery.history.get_nrows()):
        rmsds = list(refinery.history["rmsd"][i])
        rows.append(
            [str(i), str(refinery.history["num_reflections"][i])]
            + [f"{r:.5g}" for r in rmsds]
        )

    logger.info(tabulate(rows, header))
    logger.info(refinery.history.reason_for_termination)


class ScalingRefinery:
    "mixin class to add extra return method"

    def __init__(self, scaler, target, prediction_parameterisation, *args, **kwargs):
        self._target = target
        self._scaler = scaler
        self._rmsd_tolerance = scaler.params.scaling_refinery.rmsd_tolerance
        self._parameters = prediction_parameterisation

    def print_step_table(self):
        print_step_table(self)

    @property
    def rmsd_tolerance(self):
        return self._rmsd_tolerance

    def set_tolerance(self, tolerance):
        """Set a tolerance."""
        self._rmsd_tolerance = tolerance

    def test_rmsd_convergence(self):
        """Test for convergence of RMSDs"""

        # http://en.wikipedia.org/wiki/
        # Non-linear_least_squares#Convergence_criteria
        try:
            r1 = self.history["rmsd"][-1]
            r2 = self.history["rmsd"][-2]
        except IndexError:
            return False

        if r2[0] > 0:
            return abs((r2[0] - r1[0]) / r2[0]) < self._rmsd_tolerance
        else:
            return True

    def prepare_for_step(self):
        """Update the parameterisation and prepare the target function. Overwrites
        the prepare_for_step method from refinery to direct the updating away from
        the target function to the update_for_minimisation method."""

        x = self.x

        # set current parameter values
        self._parameters.set_param_vals(x)

        return

    def update_journal(self):
        """Append latest step information to the journal attributes"""

        # add step quantities to journal
        self.history.add_row()
        self.history.set_last_cell("num_reflections", self._scaler.Ih_table.size)
        self.history.set_last_cell(
            "rmsd", self._target.rmsds(self._scaler.Ih_table, self._parameters)
        )
        self.history.set_last_cell(
            "parameter_vector", self._parameters.get_param_vals()
        )
        self.history.set_last_cell("objective", self._f)
        if "gradient" in self.history:
            self.history.set_last_cell("gradient", self._g)
        return


class ScalingSimpleLBFGS(ScalingRefinery, SimpleLBFGS):
    """Adapt Refinery for L-BFGS minimiser"""

    def __init__(self, scaler, *args, **kwargs):
        logger.info("Performing a round of scaling with an LBFGS minimizer. \n")
        ScalingRefinery.__init__(self, scaler, *args, **kwargs)
        SimpleLBFGS.__init__(self, *args, **kwargs)

    def compute_functional_gradients_and_curvatures(self):
        """overwrite method to avoid calls to 'blocks' methods of target"""
        self.prepare_for_step()

        work_blocks = self._scaler.get_blocks_for_minimisation()

        f = []
        gi = []
        for block_id, block in enumerate(work_blocks):
            self._scaler.update_for_minimisation(self._parameters, block_id)
            fb, gb = self._parameters.compute_functional_gradients(block)
            f.append(fb)
            gi.append(gb)
        """task_results = easy_mp.parallel_map(
            func=self._target.compute_functional_gradients,
            iterable=blocks,
            processes=self._scaler.params.scaling_options.nproc,
            method="multiprocessing",
            preserve_exception_message=True
        )
        f, gi = zip(*task_results)"""

        f = sum(f)
        g = gi[0]
        for i in range(1, len(gi)):
            g += gi[i]

        restraints = self._parameters.compute_restraints_functional_gradients(
            self._parameters
        )

        if restraints:
            f += restraints[0]
            g += restraints[1]
        logger.debug("Functional : %s", f)
        logger.debug("Gradients : %s", list(g))
        log_memory_usage()
        logger.debug("\n")
        return f, g, None


class ScalingLstbxBuildUpMixin(ScalingRefinery):
    """Mixin class to overwrite the build_up method in AdaptLstbx"""

    def build_up(self, objective_only=False):
        "overwrite method from Adaptlstbx"
        # set current parameter values
        self.prepare_for_step()

        # Reset the state to construction time, i.e. no equations accumulated
        self.reset()

        work_blocks = self._scaler.get_blocks_for_minimisation()

        # observation terms
        if objective_only:
            for block_id, block in enumerate(work_blocks):
                self._scaler.update_for_minimisation(self._parameters, block_id)
                residuals, weights = self._parameters.compute_residuals(block)
                self.add_residuals(residuals, weights)
        else:
            self._jacobian = None

            for block_id, block in enumerate(work_blocks):
                self._scaler.update_for_minimisation(self._parameters, block_id)
                (
                    residuals,
                    jacobian,
                    weights,
                ) = self._parameters.compute_residuals_and_gradients(block)
                self.add_equations(residuals, jacobian, weights)
            """task_results = easy_mp.pool_map(
                fixed_func=self._target.compute_residuals_and_gradients,
                iterable=blocks,
                processes=self._scaler.params.scaling_options.nproc
            )
            for result in task_results:
                self.add_equations(result[0], result[1], result[2])"""

        restraints = self._parameters.compute_restraints_residuals_and_gradients(
            self._parameters
        )
        if restraints:
            if objective_only:
                self.add_residuals(restraints[0], restraints[2])
            else:
                self.add_equations(restraints[0], restraints[1], restraints[2])
        logger.debug("added equations for all blocks")
        log_memory_usage()
        logger.debug("\n")
        return


class ScalingGaussNewtonIterations(ScalingLstbxBuildUpMixin, GaussNewtonIterations):
    """Refinery implementation, using lstbx Gauss Newton iterations"""

    # defaults that may be overridden
    gradient_threshold = 1.0e-10
    step_threshold = None
    damping_value = 0.0007
    max_shift_over_esd = 15
    convergence_as_shift_over_esd = 1e-5

    def __init__(
        self,
        scaler,
        target,
        prediction_parameterisation,
        constraints_manager=None,
        log=None,
        tracking=None,
        max_iterations=20,
    ):
        logger.info("Performing a round of scaling with a Gauss-Newton minimizer.\n")
        ScalingLstbxBuildUpMixin.__init__(
            self, scaler, target, prediction_parameterisation
        )
        GaussNewtonIterations.__init__(
            self,
            target,
            prediction_parameterisation,
            constraints_manager,
            log=log,
            tracking=tracking,
            max_iterations=max_iterations,
        )


class ScalingLevenbergMarquardtIterations(
    ScalingLstbxBuildUpMixin, LevenbergMarquardtIterations
):
    """Refinery implementation, employing lstbx Levenberg Marquadt
    iterations"""

    def __init__(
        self,
        scaler,
        target,
        prediction_parameterisation,
        constraints_manager=None,
        log=None,
        tracking=None,
        max_iterations=20,
    ):
        logger.info(
            "Performing a round of scaling with a Levenberg-Marquardt minimizer.\n"
        )
        ScalingLstbxBuildUpMixin.__init__(
            self, scaler, target, prediction_parameterisation
        )
        LevenbergMarquardtIterations.__init__(
            self,
            target,
            prediction_parameterisation,
            constraints_manager,
            log=log,
            tracking=tracking,
            max_iterations=max_iterations,
        )
