"""
Refinement engine and functions for error model refinement.
"""
from __future__ import annotations

import logging

from dials.algorithms.refinement.engine import SimpleLBFGS
from dials.algorithms.scaling.error_model.error_model import (
    ErrorModelA_APM,
    ErrorModelB_APM,
    ErrorModelRegressionAPM,
)
from dials.algorithms.scaling.error_model.error_model_target import (
    ErrorModelTargetA,
    ErrorModelTargetB,
    ErrorModelTargetRegression,
)
from dials.algorithms.scaling.scaling_refiner import print_step_table

logger = logging.getLogger("dials")


def run_error_model_refinement(model, Ih_table):
    """
    Refine an error model for the input data, returning the model.

    Raises:
        ValueError: if insufficient reflections or bad refined value.
        RuntimeError: can be raised in LBFGS minimiser.
    """
    assert Ih_table.n_work_blocks == 1
    model.configure_for_refinement(Ih_table.blocked_data_list[0])
    if not model.active_parameters:
        logger.info("All error model parameters fixed, skipping refinement")
    else:
        logger.info("Performing a round of error model refinement.")
        refinery = error_model_refinery(
            model=model,
            active_parameters=model.active_parameters,
            error_model_scope=model.params,
            max_iterations=100,
        )
        if refinery:
            refinery.run()
            refinery.print_step_table()
        logger.info(model)
    model.finalise()
    return model


def error_model_refinery(model, active_parameters, error_model_scope, max_iterations):
    """Return the correct engine based on phil parameters.

    Note that here the target also takes the role of the predication
    parameterisation by implementing the set_param_vals and get_param_vals
    methods (the code is organised in this way to allow the use of the
    dials.refinement engines)."""
    if model.id_ == "basic":
        if error_model_scope.minimisation == "individual":
            return ErrorModelRefinery(
                model=model,
                max_iterations=max_iterations,
                parameters_to_refine=active_parameters,
            )
        if error_model_scope.minimisation == "regression":
            parameterisation = ErrorModelRegressionAPM(model, active_parameters)
            target = ErrorModelTargetRegression(model)
            return ErrorModelRegressionRefiner(
                model=model,
                target=target,
                prediction_parameterisation=parameterisation,
                max_iterations=max_iterations,
            )
        else:
            return None


class ErrorModelRegressionRefiner(SimpleLBFGS):

    """Use LBFGS for convenience, actually is a linear regression.

    Therefore target.predict step is unnecessary."""

    def __init__(self, model, *args, **kwargs):
        self.model = model
        self.parameterisation = kwargs["prediction_parameterisation"]
        SimpleLBFGS.__init__(self, *args, **kwargs)

    def print_step_table(self):
        print_step_table(self)

    def update_journal(self):
        """Append latest step information to the journal attributes"""

        # add step quantities to journal
        self.history.add_row()
        self.history.set_last_cell("num_reflections", self.model.n_refl)
        self.history.set_last_cell("rmsd", self._target.rmsds(self.parameterisation))
        self.history.set_last_cell(
            "parameter_vector", self.parameterisation.get_param_vals()
        )
        self.history.set_last_cell("objective", self._f)
        if "gradient" in self.history:
            self.history.set_last_cell("gradient", self._g)
        return

    def run(self):
        super().run()
        self.parameterisation.resolve_model_parameters()

    def prepare_for_step(self):
        """Update the parameterisation and prepare the target function. Overwrites
        the prepare_for_step method from refinery to direct the updating away from
        the target function to the update_for_minimisation method."""

        x = self.x

        # set current parameter values
        self.parameterisation.set_param_vals(x)

        return

    def compute_functional_gradients_and_curvatures(self):
        """overwrite method to avoid calls to 'blocks' methods of target"""
        logger.debug("Current parameters %s", [f"{i:.6f}" for i in self.x])
        self.prepare_for_step()
        self._target.predict(
            self.parameterisation
        )  # null operation as linear regression

        f, g = self._target.compute_functional_gradients(self.parameterisation)

        # restraints terms
        restraints = self._target.compute_restraints_functional_gradients(
            self.parameterisation
        )

        if restraints:
            f += restraints[0]
            g += restraints[1]
        logger.debug("Current functional %s", f)
        return f, g, None


class ErrorModelRefinery:

    """Refiner for the basic error model."""

    def __init__(self, model, parameters_to_refine, *args, **kwargs):
        self.model = model
        self.parameters_to_refine = parameters_to_refine
        self.args = args
        self.kwargs = kwargs
        self.avals = []
        self.bvals = []
        self._avals_tolerance = 0.01
        self.converged = False

    def print_step_table(self):
        pass

    def test_value_convergence(self):
        """Test for convergence of RMSDs"""

        # http://en.wikipedia.org/wiki/
        # Non-linear_least_squares#Convergence_criteria
        try:
            r1 = self.avals[-1]
            r2 = self.avals[-2]
        except IndexError:
            return False

        if r2 > 0:
            return abs((r2 - r1) / r2) < self._avals_tolerance
        else:
            return True

    def _refine_a(self):
        parameterisation = ErrorModelA_APM(self.model)
        self._refine_component(
            self.model, ErrorModelTargetA(self.model), parameterisation
        )

    def _refine_b(self):
        parameterisation = ErrorModelB_APM(self.model)
        self._refine_component(
            self.model, ErrorModelTargetB(self.model), parameterisation
        )

    def _refine_component(self, model, target, parameterisation):
        refiner = ErrorModelRegressionRefiner(
            model=model,
            target=target,
            prediction_parameterisation=parameterisation,
            *self.args,
            **self.kwargs,
        )
        refiner.run()

    def run(self):
        """Refine the model."""
        if self.parameters_to_refine == ["a", "b"]:
            for n in range(20):  # usually converges in around 5 cycles
                self._refine_a()
                # now update in model
                self.avals.append(self.model.components["a"].parameters[0])
                logger.debug(
                    "Error model refinement cycle %s: a = %s", n, self.avals[-1]
                )
                self.model.update(self.model.parameters)
                self._refine_b()
                self.bvals.append(self.model.components["b"].parameters[0])
                logger.debug(
                    "Error model refinement cycle %s: b = %s", n, self.bvals[-1]
                )
                self.model.update(self.model.parameters)
                if self.test_value_convergence():
                    self.converged = True
                    break
            if self.converged and self.avals[-1] > 0.4:
                self._refine_a()
                self.avals.append(self.model.components["a"].parameters[0])
                self.model.update(self.model.parameters)
            else:
                logger.info(
                    """
Two-parameter Error model refinement failed.
Performing error model refinement with fixed a=1.0
"""
                )
                self.model.parameters = [1.0, 0.02]
                self._refine_b()
                self.model.update(self.model.parameters)
        elif self.parameters_to_refine == ["a"]:
            self._refine_a()
            self.model.update(self.model.parameters)
        elif self.parameters_to_refine == ["b"]:
            self._refine_b()
            self.model.update(self.model.parameters)
