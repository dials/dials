""" Classes for scaling refinement engines.

Classes are inherited from the dials.refinement engine with a few
methods overwritten to use them with scaling code."""

from __future__ import absolute_import, division, print_function
import logging
from dials.algorithms.refinement.engine import (
    SimpleLBFGS,
    GaussNewtonIterations,
    LevenbergMarquardtIterations,
    LBFGScurvs,
)
from dials.algorithms.scaling.scaling_utilities import log_memory_usage
from libtbx.phil import parse
from libtbx.table_utils import simple_table
from scitbx.array_family import flex
from iotbx import merging_statistics

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
    elif engine == "LBFGScurvs":
        return ScalingLBFGScurvs(
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


def error_model_refinery(engine, target, max_iterations):
    """Return the correct engine based on phil parameters.

    Note that here the target also takes the role of the predication
    parameterisation by implementing the set_param_vals and get_param_vals
    methods (the code is organised in this way to allow the use of the
    dials.refinement engines)."""
    if engine == "SimpleLBFGS":
        return ErrorModelSimpleLBFGS(
            target=target,
            prediction_parameterisation=target,
            max_iterations=max_iterations,
        )
    """elif engine == 'LBFGScurvs':
    assert 0, 'LBFGS with curvatures not yet implemented'
  elif engine == 'GaussNewton':
    return ScalingGaussNewtonIterations(target=target,
      prediction_parameterisation=prediction_parameterisation,
      max_iterations=max_iterations)
  elif engine == 'LevMar':
    return ScalingLevenbergMarquardtIterations(target=target,
      prediction_parameterisation=prediction_parameterisation,
      max_iterations=max_iterations)"""


def print_step_table(refinery):
    """print useful output about refinement steps in the form of a simple table"""

    logger.info("\nRefinement steps:")

    header = ["Step", "Nref"]
    for (name, units) in zip(refinery._target.rmsd_names, refinery._target.rmsd_units):
        header.append(name + "\n(" + units + ")")

    rows = []
    for i in range(refinery.history.get_nrows()):
        rmsds = [r for r in refinery.history["rmsd"][i]]
        rows.append(
            [str(i), str(refinery.history["num_reflections"][i])]
            + ["%.5g" % r for r in rmsds]
        )

    st = simple_table(rows, header)
    logger.info(st.format())
    logger.info(refinery.history.reason_for_termination)


class ErrorModelRefinery(object):
    """Mixin class to add extra return method."""

    def __init__(self, error_model_target):
        self.error_model_target = error_model_target

    def return_error_model(self):
        """Set error manager parameters and return error manager."""
        print_step_table(self)
        self._target.error_model.refined_parameters = self._target.x
        # logger.info("\nMinimised error model with parameters {0:.5f} and {1:.5f}. {sep}"
        #      .format(self._target.x[0], abs(self._target.x[1]), sep='\n'))
        self._target.error_model.intensities = (
            self._target.error_model.Ih_table.intensities
        )
        self._target.error_model.sigmaprime = self._target.error_model.sigmaprime
        self._target.error_model.inverse_scale_factors = (
            self._target.error_model.Ih_table.inverse_scale_factors
        )
        self._target.error_model.clear_Ih_table()
        return self._target.error_model


class ScalingRefinery(object):
    "mixin class to add extra return method"

    def __init__(self, scaler, target, prediction_parameterisation, *args, **kwargs):
        self._target = target
        self._scaler = scaler
        self._rmsd_tolerance = scaler.params.scaling_refinery.rmsd_tolerance
        self._parameters = prediction_parameterisation

    def test_rmsd_convergence(self):
        """Test for convergence of RMSDs"""

        # http://en.wikipedia.org/wiki/
        # Non-linear_least_squares#Convergence_criteria
        try:
            r1 = self.history["rmsd"][-1]
            r2 = self.history["rmsd"][-2]
        except IndexError:
            return False

        tests = [
            abs((r2[0] - r1[0]) / r2[0]) < self._rmsd_tolerance if r2[0] > 0 else True
        ]

        return all(tests)

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

    def return_scaler(self):
        """return scaler method"""
        from dials.algorithms.scaling.scaler import MultiScalerBase

        print_step_table(self)

        if self._scaler.id_ == "single":
            if self._parameters.apm_list[0].var_cov_matrix:
                self._scaler.update_var_cov(self._parameters.apm_list[0])
                self._scaler.experiment.scaling_model.set_scaling_model_as_scaled()
        elif self._scaler.id_ == "multi" or self._scaler.id_ == "target":
            if self._parameters.apm_list[0].var_cov_matrix:  # test if has been set
                for i, scaler in enumerate(self._scaler.active_scalers):
                    scaler.update_var_cov(self._parameters.apm_list[i])
                    scaler.experiment.scaling_model.set_scaling_model_as_scaled()

        if not isinstance(self._scaler, MultiScalerBase):
            self._scaler.experiment.scaling_model.normalise_components()

        logger.debug("\n" + str(self._scaler.experiment.scaling_model))

        if self._scaler.Ih_table.free_Ih_table:
            i_obs = self._scaler.Ih_table.as_miller_array(
                self._scaler.experiment.crystal.get_unit_cell(),
                return_free_set_data=True,
            )
            res = merging_statistics.dataset_statistics(
                i_obs=i_obs,
                n_bins=20,
                anomalous=False,
                sigma_filtering=None,
                use_internal_variance=False,
                eliminate_sys_absent=False,
                cc_one_half_method="sigma_tau",
            )
            free_rpim = res.overall.r_pim
            free_cc12 = res.overall.cc_one_half
            ccs = flex.double([b.cc_one_half for b in res.bins])
            n_refl = flex.double([b.n_obs for b in res.bins])
            n_tot = sum(n_refl)
            free_wcc12 = sum(ccs * n_refl / n_tot)

            i_obs = self._scaler.Ih_table.as_miller_array(
                self._scaler.experiment.crystal.get_unit_cell()
            )
            res = merging_statistics.dataset_statistics(
                i_obs=i_obs,
                n_bins=20,
                anomalous=False,
                sigma_filtering=None,
                use_internal_variance=False,
                eliminate_sys_absent=False,
                cc_one_half_method="sigma_tau",
            )
            work_rpim = res.overall.r_pim
            work_cc12 = res.overall.cc_one_half
            ccs = flex.double([b.cc_one_half for b in res.bins])
            n_refl = flex.double([b.n_obs for b in res.bins])
            n_tot = sum(n_refl)
            work_wcc12 = sum(ccs * n_refl / n_tot)

            rpim_gap = free_rpim - work_rpim
            cc12_gap = work_cc12 - free_cc12
            wcc12_gap = work_wcc12 - free_wcc12
            self._scaler.final_rmsds = [
                work_rpim,
                free_rpim,
                rpim_gap,
                work_cc12,
                free_cc12,
                cc12_gap,
                work_wcc12,
                free_wcc12,
                wcc12_gap,
            ]
            header = ["", "Work", "Free", "Gap"]
            rows = [
                [
                    "Rpim",
                    str(round(work_rpim, 5)),
                    str(round(free_rpim, 5)),
                    str(round(rpim_gap, 5)),
                ],
                [
                    "CC1/2",
                    str(round(work_cc12, 5)),
                    str(round(free_cc12, 5)),
                    str(round(cc12_gap, 5)),
                ],
                [
                    "CC1/2 (weighted-avg)",
                    str(round(work_wcc12, 5)),
                    str(round(free_wcc12, 5)),
                    str(round(wcc12_gap, 5)),
                ],
            ]
            logger.info(
                """\nWork/Free set quality indicators:
(CC1/2 calculated using the sigma-tau method, weighted-avg is the
average CC1/2 over resolution bins, weighted by n_obs per bin.)
Gaps are defined as Rfree-Rwork and CCWork-CCfree."""
            )
            st = simple_table(rows, header)
            logger.info(st.format())


class ScalingSimpleLBFGS(ScalingRefinery, SimpleLBFGS):
    """Adapt Refinery for L-BFGS minimiser"""

    def __init__(self, scaler, *args, **kwargs):
        logger.info("Performing a round of scaling with an LBFGS minimizer. \n")
        ScalingRefinery.__init__(self, scaler, *args, **kwargs)
        SimpleLBFGS.__init__(self, *args, **kwargs)

    def compute_functional_gradients_and_curvatures(self):
        """overwrite method to avoid calls to 'blocks' methods of target"""
        self.prepare_for_step()

        if self._scaler.Ih_table.free_Ih_table:
            blocks = self._scaler.Ih_table.blocked_data_list[:-1]
            free_block_id = len(self._scaler.Ih_table.blocked_data_list) - 1
            self._scaler.update_for_minimisation(self._parameters, free_block_id)
            self._scaler.clear_memory_from_derivs(free_block_id)
        else:
            blocks = self._scaler.Ih_table.blocked_data_list

        # if self._scaler.params.scaling_options.nproc > 1:
        f = []
        gi = []
        for block_id, block in enumerate(blocks):
            self._scaler.update_for_minimisation(self._parameters, block_id)
            fb, gb = self._target.compute_functional_gradients(block)
            f.append(fb)
            gi.append(gb)
            self._scaler.clear_memory_from_derivs(block_id)
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
        """else:
      f = 0.0
      g = None
      for block in blocks:
        fi, gi = self._target.compute_functional_gradients(block)
        f += fi
        if g:
          g += gi
        else:
          g = gi"""

        restraints = self._target.compute_restraints_functional_gradients_and_curvatures(
            self._parameters
        )

        if restraints:
            f += restraints[0]
            g += restraints[1]
        logger.debug("Functional : %s" % f)
        logger.debug("Gradients : %s" % list(g))
        log_memory_usage()
        logger.debug("\n")
        return f, g, None


class ScalingLBFGScurvs(ScalingRefinery, LBFGScurvs):
    """Adapt Refinery for L-BFGS minimiser"""

    def __init__(self, scaler, *args, **kwargs):
        logger.info("Performing a round of scaling with an LBFGS minimizer. \n")
        LBFGScurvs.__init__(self, *args, **kwargs)
        ScalingRefinery.__init__(self, scaler, *args, **kwargs)
        self._target.curvatures = True

    def compute_functional_gradients_and_curvatures(self):
        """overwrite method to avoid calls to 'blocks' methods of target"""

        self.prepare_for_step()

        f, g, curv = self._target.compute_functional_gradients_and_curvatures()

        # restraints terms
        restraints = (
            self._target.compute_restraints_functional_gradients_and_curvatures()
        )

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
        logger.debug("Current parameters %s" % ["%.6f" % i for i in self.x])
        self.prepare_for_step()

        f, g, _ = self._target.compute_functional_gradients_and_curvatures()

        # restraints terms
        restraints = (
            self._target.compute_restraints_functional_gradients_and_curvatures()
        )

        if restraints:
            f += restraints[0]
            g += restraints[1]
        logger.debug("Current functional %s" % f)
        return f, g, None


class ScalingLstbxBuildUpMixin(ScalingRefinery):
    """Mixin class to overwrite the build_up method in AdaptLstbx"""

    def build_up(self, objective_only=False):
        "overwrite method from Adaptlstbx"
        # set current parameter values
        self.prepare_for_step()

        # Reset the state to construction time, i.e. no equations accumulated
        self.reset()

        if self._scaler.Ih_table.free_Ih_table:
            blocks = self._scaler.Ih_table.blocked_data_list[:-1]
            free_block_id = len(self._scaler.Ih_table.blocked_data_list) - 1
            self._scaler.update_for_minimisation(self._parameters, free_block_id)
            self._scaler.clear_memory_from_derivs(free_block_id)
        else:
            blocks = self._scaler.Ih_table.blocked_data_list

        # observation terms
        if objective_only:
            # if self._scaler.params.scaling_options.nproc > 1: #no mp option yet
            for block_id, block in enumerate(blocks):
                self._scaler.update_for_minimisation(self._parameters, block_id)
                residuals, weights = self._target.compute_residuals(block)
                self.add_residuals(residuals, weights)
        else:
            # if self._scaler.params.scaling_options.nproc: #no mp option yet

            self._jacobian = None

            for block_id, block in enumerate(blocks):
                self._scaler.update_for_minimisation(self._parameters, block_id)
                residuals, jacobian, weights = self._target.compute_residuals_and_gradients(
                    block
                )
                self.add_equations(residuals, jacobian, weights)
            """task_results = easy_mp.pool_map(
        fixed_func=self._target.compute_residuals_and_gradients,
        iterable=blocks,
        processes=self._scaler.params.scaling_options.nproc
        )
      for result in task_results:
        self.add_equations(result[0], result[1], result[2])"""

        restraints = self._target.compute_restraints_residuals_and_gradients(
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
