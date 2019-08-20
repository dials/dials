#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""Contains classes for refinement engines. Refinery is the shared interface,
LevenbergMarquardtIterations, GaussNewtonIterations, SimpleLBFGS and LBFGScurvs
are the current concrete implementations"""

from __future__ import absolute_import, division, print_function

import copy
import logging
import json

import libtbx
from dials.algorithms.refinement import DialsRefineRuntimeError
from libtbx import easy_mp
from libtbx.phil import parse
from scitbx import lbfgs
from scitbx.array_family import flex
from scitbx.lstbx import normal_eqns, normal_eqns_solving

logger = logging.getLogger(__name__)


# termination reason strings
TARGET_ACHIEVED = "RMSD target achieved"
RMSD_CONVERGED = "RMSD no longer decreasing"
STEP_TOO_SMALL = "Step too small"
OBJECTIVE_INCREASE = "Refinement failure: objective increased"
MAX_ITERATIONS = "Reached maximum number of iterations"
MAX_TRIAL_ITERATIONS = "Reached maximum number of consecutive unsuccessful trial steps"
DOF_TOO_LOW = "Not enough degrees of freedom to refine"

refinery_phil_str = """
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
"""
refinery_phil_scope = parse(refinery_phil_str)


class Journal(dict):
    """Container in which to store information about refinement history.

    This is simply a dict but provides some extra methods for access that
    maintain values as columns in a table. Refinery classes will use these methods
    while entering data to ensure the table remains consistent. Methods inherited
    from dict are not hidden for ease of use of this object when returned to the
    user."""

    reason_for_termination = None
    _nrows = 0

    def get_nrows(self):
        return self._nrows

    def add_column(self, key):
        """Add a new column named by key"""
        self[key] = [None] * self._nrows

    def add_row(self):
        """Add an element to the end of each of the columns. Fail if any columns
        are the wrong length"""

        for k in self:
            assert len(self[k]) == self._nrows
            self[k].append(None)
        self._nrows += 1

    def del_last_row(self):
        """Delete the last element from the each of the columns. Fail if any columns
        are the wrong length"""

        if self._nrows == 0:
            return None
        for k in self:
            assert len(self[k]) == self._nrows
            self[k].pop()
        self._nrows -= 1

    def set_last_cell(self, key, value):
        """Set last cell in column given by key to value. Fail if the column is the
        wrong length"""

        assert len(self[key]) == self._nrows
        self[key][-1] = value

    def to_json_file(self, filename):
        d = {"attributes": self.__dict__, "data": dict(self)}
        with open(filename, "w") as f:
            json.dump(d, f)

    @classmethod
    def from_json_file(cls, filename):
        with open(filename, "r") as f:
            d = json.load(f)
        j = cls()
        j.update(d["data"])
        for key in d["attributes"]:
            setattr(j, key, d["attributes"][key])
        return j


class Refinery(object):
    """Interface for Refinery objects. This should be subclassed and the run
    method implemented."""

    # NOTES. A Refinery is initialised with a Target function. The target
    # function already contains a ReflectionManager (which holds the data) so
    # there's no need to pass the data in here. In fact the Target
    # class does the bulk of the work, as it also does the reflection prediction
    # to get the updated predictions on each cycle. This should make some sense
    # as the target function is inextricably linked to the space in which
    # predictions are made (e.g. detector space, phi), so it is not general
    # enough to sit abstractly above the prediction.

    # This keeps the Refinery simple and able to be focused only on generic
    # features of managing a refinement run, like reporting results and checking
    # termination criteria.

    # The prediction values come from a PredictionParameterisation object.
    # This is also referred to by the Target function, but it makes sense for
    # Refinery to be able to refer to it directly. So refinery should keep a
    # separate link to its PredictionParameterisation.

    def __init__(
        self,
        target,
        prediction_parameterisation,
        constraints_manager=None,
        log=None,
        tracking=None,
        max_iterations=None,
    ):

        # reference to PredictionParameterisation, Target and ConstraintsManager
        # objects
        self._parameters = prediction_parameterisation
        self._target = target
        self._constr_manager = constraints_manager

        # initial parameter values
        self.x = flex.double(self._parameters.get_param_vals())
        if self._constr_manager is not None:
            self.x = self._constr_manager.constrain_parameters(self.x)
        self.old_x = None

        # undefined initial functional and gradients values
        self._f = None
        self._g = None
        self._jacobian = None

        # filename for an optional log file
        self._log = log

        self._target_achieved = False
        self._max_iterations = max_iterations

        # attributes for journalling functionality, based on lstbx's
        # journaled_non_linear_ls class
        if tracking is None:
            # set default tracking
            tracking = refinery_phil_scope.extract().refinery.journal
        self.history = Journal()
        self.history.add_column("num_reflections")
        self.history.add_column("objective")  # flex.double()
        if tracking.track_gradient:
            self.history.add_column("gradient")
        self.history.add_column("gradient_norm")  # flex.double()
        if tracking.track_parameter_correlation:
            self.history.add_column("parameter_correlation")
        if tracking.track_step:
            self.history.add_column("solution")
        if tracking.track_out_of_sample_rmsd:
            self.history.add_column("out_of_sample_rmsd")
        self.history.add_column("solution_norm")  # flex.double()
        self.history.add_column("parameter_vector")
        self.history.add_column("parameter_vector_norm")  # flex.double()
        self.history.add_column("rmsd")
        if tracking.track_condition_number:
            self.history.add_column("condition_number")

        # number of processes to use, for engines that support multiprocessing
        self._nproc = 1

        self.prepare_for_step()

    def get_num_steps(self):
        return self.history.get_nrows() - 1

    def prepare_for_step(self):
        """Update the parameterisation and prepare the target function"""

        x = self.x
        if self._constr_manager is not None:
            x = self._constr_manager.expand_parameters(x)

        # set current parameter values
        self._parameters.set_param_vals(x)

        # do reflection prediction
        self._target.predict()

    def update_journal(self):
        """Append latest step information to the journal attributes"""

        # add step quantities to journal
        self.history.add_row()
        self.history.set_last_cell("num_reflections", self._target.get_num_matches())
        self.history.set_last_cell("rmsd", self._target.rmsds())
        self.history.set_last_cell(
            "parameter_vector", self._parameters.get_param_vals()
        )
        self.history.set_last_cell("objective", self._f)
        if "gradient" in self.history:
            self.history.set_last_cell("gradient", self._g)
        if "parameter_correlation" in self.history and self._jacobian is not None:
            resid_names = [s.replace("RMSD_", "") for s in self._target.rmsd_names]
            # split Jacobian into dense matrix blocks corresponding to each residual
            jblocks = self.split_jacobian_into_blocks()
            corrmats = {}
            for r, j in zip(resid_names, jblocks):
                corrmats[r] = self._packed_corr_mat(j)
            self.history.set_last_cell("parameter_correlation", corrmats)
        if "condition_number" in self.history and self._jacobian is not None:
            self.history.set_last_cell(
                "condition_number", self.jacobian_condition_number()
            )
        if "out_of_sample_rmsd" in self.history:
            preds = self._target.predict_for_free_reflections()
            self.history.set_last_cell(
                "out_of_sample_rmsd", self._target.rmsds_for_reflection_table(preds)
            )

    def split_jacobian_into_blocks(self):
        """Split the Jacobian into blocks each corresponding to a separate
        residual"""

        nblocks = len(self._target.rmsd_names)

        try:
            # The Jacobian might be a sparse matrix
            j = self._jacobian.as_dense_matrix()
        except AttributeError:
            j = self._jacobian

        nr, nc = j.all()
        nr_block = int(nr / nblocks)
        row_start = [e * nr_block for e in range(nblocks)]
        blocks = [j.matrix_copy_block(rs, 0, nr_block, nc) for rs in row_start]

        return blocks

    @staticmethod
    def _packed_corr_mat(m):
        """Return a 1D flex array containing the upper diagonal values of the
        correlation matrix calculated between columns of 2D matrix m"""

        nr, nc = m.all()

        try:  # convert a flex.double matrix to sparse
            from scitbx import sparse

            m2 = sparse.matrix(nr, nc)
            m2.assign_block(m, 0, 0)
            m = m2
        except AttributeError:
            pass  # assume m is already scitbx_sparse_ext.matrix

        packed_len = (m.n_cols * (m.n_cols + 1)) // 2
        i = 0
        tmp = flex.double(packed_len)
        for col1 in range(m.n_cols):
            for col2 in range(col1, m.n_cols):
                if col1 == col2:
                    tmp[i] = 1.0
                else:
                    # Avoid spuriously high correlation between a column that should be
                    # zero (such as the gradient of X residuals wrt the Shift2 parameter)
                    # and another column (such as the gradient of X residuals wrt the
                    # Dist parameter) by rounding values to 15 places. It seems that such
                    # spurious correlations may occur in cases where gradients are
                    # calculated to be zero by matrix operations, rather than set to zero.
                    v1 = m.col(col1).as_dense_vector().round(15)
                    v2 = m.col(col2).as_dense_vector().round(15)
                    tmp[i] = flex.linear_correlation(v1, v2).coefficient()
                i += 1

        return tmp

    def get_correlation_matrix_for_step(self, step):
        """For each type of residual (e.g. X, Y, Phi), decompress and return the
        full 2D correlation matrix between columns of the Jacobian that was
        stored in the journal at the given step number. If not available, return
        None"""

        if "parameter_correlation" not in self.history:
            return None
        try:
            packed_mats = self.history["parameter_correlation"][step]
        except IndexError:
            return None
        if packed_mats is None:
            return None

        packed_mats = copy.deepcopy(packed_mats)

        nparam = len(self._parameters)

        for k, v in packed_mats.items():
            corr_mat = flex.double(flex.grid(nparam, nparam))
            i = 0
            for row in range(nparam):
                for col in range(row, nparam):
                    corr_mat[row, col] = v[i]
                    i += 1
            corr_mat.matrix_copy_upper_to_lower_triangle_in_place()
            packed_mats[k] = corr_mat
        return packed_mats

    def jacobian_condition_number(self):
        """Calculate the condition number of the Jacobian, for tracking in the
        refinement journal, if requested. The condition number of a matrix A is
        defined as cond(A) = ||A|| ||inv(A)||. For a rectangular matrix the inverse
        operation refers to the Moore-Penrose pseudoinverse. Various matrix norms
        can be used, resulting in numerically different condition numbers, however
        the 2-norm is commonly used. In that case, the definition is equivalent
        to the ratio of the largest to smallest singular values of the matrix:
        cond(A) = sig_(A) / sig_min(A). That is the calculation that is performed
        here.

        The condition number is a measure of how accurate the solution x to the
        equation Ax = b will be. Essentially it measures how errors are amplified
        through the linear equation. The condition number is large in the case that
        the columns of A are nearly linearly-dependent (and infinite for a singular
        matrix). We use it here then to detect situations where the correlation
        between effects of different parameter shifts becomes large and therefore
        refinement is problematic.

        Note, the Jacobian used here does not include any additional rows due to
        restraints terms that might be applied, or any parameter reduction due to
        constraints. Therefore this condition number relates to the pure linearised
        (Gauss-Newton) step, which might not actually be what the refinement engine
        uses. It can be indicative of issues in the fundamental set up of the least
        squares problem, even if these issues are avoided in practice (e.g. by
        use of an algorithm like Levenberg-Marquardt, inclusion of restraints or
        parameter reduction).
        """
        try:
            # The Jacobian might be a sparse matrix
            j = self._jacobian.as_dense_matrix().deep_copy()
        except AttributeError:
            j = self._jacobian.deep_copy()

        from scitbx.linalg.svd import real as svd_real

        svd = svd_real(j, False, False)

        # The condition number is the ratio of the largest to the smallest singular
        # values of the matrix
        return max(svd.sigma) / min(svd.sigma)

    def test_for_termination(self):
        """Return True if refinement should be terminated"""

        # Basic version delegate to the Target class. Derived classes may
        # implement other termination criteria
        self._target_achieved = self._target.achieved()

        return self._target_achieved

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
            abs((e[1] - e[0]) / e[1]) < 0.0001 if e[1] > 0 else True
            for e in zip(r1, r2)
        ]

        return all(tests)

    def test_objective_increasing_but_not_nref(self):
        """Test for an increase in the objective value between steps. This
        could be caused simply by the number of matches between observations
        and predictions increasing. However, if the number of matches stayed
        the same or reduced then this is a bad sign."""

        try:
            l1 = self.history["objective"][-1]
            l2 = self.history["objective"][-2]
            n1 = self.history["num_reflections"][-1]
            n2 = self.history["num_reflections"][-2]
        except IndexError:
            return False

        return l1 > l2 and n1 <= n2

    def set_nproc(self, nproc):
        """Set number of processors for multiprocessing. Override in derived classes
        if a policy dictates that this must not be user-controlled"""
        self._nproc = nproc

    def run(self):
        """
        To be implemented by derived class. It is expected that each step of
        refinement be preceeded by a call to prepare_for_step and followed by
        calls to update_journal and test_for_termination (in that order).
        """

        # Specify a minimizer and its parameters, and run
        raise NotImplementedError()


class DisableMPmixin(object):
    """A mixin class that disables setting of nproc for multiprocessing"""

    def set_nproc(self, nproc):
        if nproc != 1:
            raise NotImplementedError()


class AdaptLbfgs(Refinery):
    """Adapt Refinery for L-BFGS minimiser"""

    def __init__(self, *args, **kwargs):

        Refinery.__init__(self, *args, **kwargs)

        self._termination_params = lbfgs.termination_parameters(
            max_iterations=self._max_iterations
        )

        from six.moves import cStringIO as StringIO

        self._log_string = StringIO

    def compute_functional_and_gradients(self):
        L, dL_dp, _ = self.compute_functional_gradients_and_curvatures()
        self._f = L
        self._g = dL_dp
        return self._f, self._g

    def compute_functional_gradients_and_curvatures(self):

        self.prepare_for_step()

        # observation terms
        blocks = self._target.split_matches_into_blocks(nproc=self._nproc)
        if self._nproc > 1:
            task_results = easy_mp.parallel_map(
                func=self._target.compute_functional_gradients_and_curvatures,
                iterable=blocks,
                processes=self._nproc,
                method="multiprocessing",
                preserve_exception_message=True,
            )

        else:
            task_results = [
                self._target.compute_functional_gradients_and_curvatures(block)
                for block in blocks
            ]

        # reduce blockwise results
        flist, glist, clist = zip(*task_results)
        f = sum(flist)
        g = [sum(g) for g in zip(*glist)]
        c = [sum(c) for c in zip(*clist)]

        # restraints terms
        restraints = (
            self._target.compute_restraints_functional_gradients_and_curvatures()
        )

        if restraints:
            f += restraints[0]
            g = [a + b for a, b in zip(g, restraints[1])]
            c = [a + b for a, b in zip(c, restraints[2])]

        # compact and reorder according to the constraints
        if self._constr_manager is not None:
            g = self._constr_manager.constrain_gradient_vector(g)
            c = self._constr_manager.constrain_gradient_vector(c)

        return f, flex.double(g), flex.double(c)

    def callback_after_step(self, minimizer):
        """
        Do journalling, evaluate rmsds and return True if the target is
        reached to terminate the refinement.
        """

        self.update_journal()
        logger.debug("Step %d", self.history.get_nrows() - 1)

        if self.test_for_termination():
            self.history.reason_for_termination = TARGET_ACHIEVED
            return True

        if self.test_rmsd_convergence():
            self.history.reason_for_termination = RMSD_CONVERGED
            return True

        return False

    def run_lbfgs(self, curvatures=False):
        """
        Run the minimiser, keeping track of its log.
        """

        ref_log = self._log_string()
        if curvatures:
            self.diag_mode = "always"
        self.minimizer = lbfgs.run(
            target_evaluator=self,
            termination_params=self._termination_params,
            log=ref_log,
        )

        log = ref_log.getvalue()
        if self._log:
            with open(self._log, "a") as f:
                f.write(log)
        ref_log.close()

        pos = log.rfind("lbfgs minimizer stop: ")
        if pos >= 0:
            msg = log[pos:].splitlines()[0]
            if self.history.reason_for_termination:
                self.history.reason_for_termination += "\n"
                self.history.reason_for_termination += msg
            else:
                self.history.reason_for_termination = msg

        if self.minimizer.error:
            self.history.reason_for_termination = self.minimizer.error


class SimpleLBFGS(AdaptLbfgs):
    """Refinery implementation, using cctbx LBFGS with basic settings"""

    def run(self):
        return self.run_lbfgs(curvatures=False)


class LBFGScurvs(AdaptLbfgs):
    """Refinery implementation using cctbx LBFGS with curvatures"""

    def run(self):
        return self.run_lbfgs(curvatures=True)

    def compute_functional_gradients_diag(self):

        L, dL_dp, curvs = self.compute_functional_gradients_and_curvatures()
        self._f = L
        self._g = dL_dp

        # Curvatures of zero will cause a crash, because their inverse is taken.
        assert curvs.all_gt(0.0)

        diags = 1.0 / curvs

        msg = "  curv: " + "%.5f " * len(tuple(curvs))
        logger.debug(msg, *curvs)

        return self._f, self._g, diags


class AdaptLstbx(Refinery, normal_eqns.non_linear_ls, normal_eqns.non_linear_ls_mixin):
    """Adapt Refinery for lstbx"""

    def __init__(
        self,
        target,
        prediction_parameterisation,
        constraints_manager=None,
        log=None,
        tracking=None,
        max_iterations=None,
    ):

        Refinery.__init__(
            self,
            target,
            prediction_parameterisation,
            constraints_manager,
            log=log,
            tracking=tracking,
            max_iterations=max_iterations,
        )

        # required for restart to work (do I need that method?)
        self.x_0 = self.x.deep_copy()

        # keep attribute for the Cholesky factor required for ESD calculation
        self.cf = None

        normal_eqns.non_linear_ls.__init__(self, n_parameters=len(self.x))

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
            residuals, weights = self._target.compute_residuals()
            self.add_residuals(residuals, weights)
        else:
            blocks = self._target.split_matches_into_blocks(nproc=self._nproc)

            if self._nproc > 1:

                # ensure the jacobian is not tracked
                self._jacobian = None

                # processing functions
                def task_wrapper(block):
                    residuals, jacobian, weights = self._target.compute_residuals_and_gradients(
                        block
                    )
                    return dict(residuals=residuals, jacobian=jacobian, weights=weights)

                def callback_wrapper(result):
                    j = result["jacobian"]
                    if self._constr_manager is not None:
                        j = self._constr_manager.constrain_jacobian(j)
                    self.add_equations(result["residuals"], j, result["weights"])
                    # no longer need the result
                    result["residuals"] = None
                    result["jacobian"] = None
                    result["weights"] = None
                    return

                easy_mp.parallel_map(
                    func=task_wrapper,
                    iterable=blocks,
                    processes=self._nproc,
                    callback=callback_wrapper,
                    method="multiprocessing",
                    preserve_exception_message=True,
                )

            else:
                for block in blocks:
                    residuals, self._jacobian, weights = self._target.compute_residuals_and_gradients(
                        block
                    )
                    j = self._jacobian
                    if self._constr_manager is not None:
                        j = self._constr_manager.constrain_jacobian(j)
                    self.add_equations(residuals, j, weights)

        # restraints terms
        restraints = self._target.compute_restraints_residuals_and_gradients()
        if restraints:
            if objective_only:
                self.add_residuals(restraints[0], restraints[2])
            else:
                j = restraints[1]
                if self._constr_manager is not None:
                    j = self._constr_manager.constrain_jacobian(j)
                self.add_equations(restraints[0], j, restraints[2])

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

    def calculate_esds(self):
        """Calculate ESDs of parameters"""

        # it is possible to get here with zero steps taken by the minimiser. For
        # example by failing for the MAX_TRIAL_ITERATIONS reason before any forward
        # steps are taken with the LevMar engine. If so the below is invalid,
        # so return early
        if self.history.get_nrows() == 0:
            return None

        if self.cf is None:
            return None

        # if constraints were used then the normal matrix has fewer rows/columns
        # than the number of expanded parameters. At the moment, do not support
        # this calculation when constraints were used
        if self._constr_manager is not None:
            return None

        # invert normal matrix from N^-1 = (U^-1)(U^-1)^T
        cf_inv = self.cf.matrix_packed_u_as_upper_triangle().matrix_inversion()
        nm_inv = cf_inv.matrix_multiply_transpose(cf_inv)

        # keep the estimated parameter variance-covariance matrix
        self.parameter_var_cov = self.history["reduced_chi_squared"][-1] * nm_inv
        # send this back to the models to calculate their uncertainties
        self._parameters.calculate_model_state_uncertainties(self.parameter_var_cov)

        # send parameter variances back to the parameter classes
        # themselves, for reporting purposes and for building restraints
        # based on existing parameterisations.
        s2 = self.parameter_var_cov.matrix_diagonal()
        assert s2.all_ge(0.0)
        s = flex.sqrt(s2)
        self._parameters.set_param_esds(s)

    def _print_normal_matrix(self):
        """Print the full normal matrix at the current step. For debugging only"""
        logger.debug("The normal matrix for the current step is:")
        logger.debug(
            self.normal_matrix_packed_u()
            .matrix_packed_u_as_symmetric()
            .as_scitbx_matrix()
            .matlab_form(format=None, one_row_per_line=True)
        )
        logger.debug("\n")


class GaussNewtonIterations(AdaptLstbx, normal_eqns_solving.iterations):
    """Refinery implementation, using lstbx Gauss Newton iterations"""

    # defaults that may be overridden
    gradient_threshold = 1.0e-10
    step_threshold = None
    damping_value = 0.0007
    max_shift_over_esd = 15
    convergence_as_shift_over_esd = 1e-5

    def __init__(
        self,
        target,
        prediction_parameterisation,
        constraints_manager=None,
        log=None,
        tracking=None,
        max_iterations=20,
        **kwds
    ):

        AdaptLstbx.__init__(
            self,
            target,
            prediction_parameterisation,
            constraints_manager,
            log=log,
            tracking=tracking,
            max_iterations=max_iterations,
        )

        # add an attribute to the journal
        self.history.add_column("reduced_chi_squared")  # flex.double()

        # adopt any overrides of the defaults above
        libtbx.adopt_optional_init_args(self, kwds)

    def run(self):
        self.n_iterations = 0

        # prepare for first step
        self.build_up()

        # return early if refinement is not possible
        if self.dof < 1:
            self.history.reason_for_termination = DOF_TOO_LOW
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

            # standard journalling
            self.update_journal()
            logger.debug("Step %d", self.history.get_nrows() - 1)

            # add cached items to the journal
            self.history.set_last_cell("parameter_vector_norm", pvn)
            self.history.set_last_cell("gradient_norm", gn)

            # extra journalling post solve
            if "solution" in self.history:
                self.history.set_last_cell("solution", self.actual.step().deep_copy())
            self.history.set_last_cell("solution_norm", self.step().norm())
            self.history.set_last_cell("reduced_chi_squared", self.chi_sq())

            # test termination criteria
            if self.test_for_termination():
                self.history.reason_for_termination = TARGET_ACHIEVED
                break

            if self.test_rmsd_convergence():
                self.history.reason_for_termination = RMSD_CONVERGED
                break

            if self.had_too_small_a_step():
                self.history.reason_for_termination = STEP_TOO_SMALL
                break

            if self.test_objective_increasing_but_not_nref():
                self.history.reason_for_termination = OBJECTIVE_INCREASE
                if self.step_backward():
                    self.history.reason_for_termination += (
                        ". Parameters set back one step"
                    )
                self.prepare_for_step()
                break

            if self.n_iterations == self._max_iterations:
                self.history.reason_for_termination = MAX_ITERATIONS
                break

            # prepare for next step
            self.step_forward()
            self.n_iterations += 1
            self.build_up()

        self.set_cholesky_factor()
        self.calculate_esds()


class LevenbergMarquardtIterations(GaussNewtonIterations):
    """Refinery implementation, employing lstbx Levenberg Marquadt
    iterations"""

    tau = 1e-3

    @property
    def mu(self):
        return self._mu

    @mu.setter
    def mu(self, value):
        self._mu = value

    def setup_mu(self):
        """Setup initial value for mu"""
        a = self.normal_matrix_packed_u()
        self.mu = self.tau * flex.max(a.matrix_packed_u_diagonal())

    def add_constant_to_diagonal(self, mu):
        """Add the constant value mu to the diagonal of the normal matrix"""
        a = self.normal_matrix_packed_u()
        a.matrix_packed_u_diagonal_add_in_place(self.mu)

    def report_progress(self, objective):
        """Callback within the refinement main loop that can be overridden to
        report the value of the objective function (and possibly) other details for
        long-running methods"""
        pass

    def _run_core(self):

        # add an attribute to the journal
        self.history.add_column("mu")
        self.history.add_column("nu")

        # FIXME need a much neater way of doing this stuff through
        # inheritance
        # set max iterations if not already.
        if self._max_iterations is None:
            self._max_iterations = 100

        self.n_iterations = 0
        nu = 2
        self.build_up()

        # early test for linear independence, require all right hand side elements to be non-zero
        RHS = self.step_equations().right_hand_side()
        if RHS.count(0.0) > 0:
            raise DialsRefineRuntimeError(
                r"""There is at least one normal equation with a right hand side of zero,
meaning that the parameters are not all independent, and there is no unique
solution.  Mathematically, some kind of row reduction needs to be performed
before this can be solved."""
            )

        # return early if refinement is not possible
        if self.dof < 1:
            self.history.reason_for_termination = DOF_TOO_LOW
            return

        self.setup_mu()

        while True:

            # set functional and gradients for the step
            self._f = self.objective()
            self._g = -self.opposite_of_gradient()

            # cache some items for the journal prior to solve
            pvn = self.parameter_vector_norm()
            gn = self.opposite_of_gradient().norm_inf()

            self.add_constant_to_diagonal(self.mu)

            # solve the normal equations
            self.solve()

            # keep the cholesky factor for ESD calculation if we end this step. Doing
            # it here ensures the normal equations are solved (cholesky_factor_packed_u
            # can only be called if that is the case)
            self.set_cholesky_factor()

            # standard journalling
            self.update_journal()
            logger.debug("Step %d", self.history.get_nrows() - 1)

            # add cached items to the journal
            self.history.set_last_cell("parameter_vector_norm", pvn)
            self.history.set_last_cell("gradient_norm", gn)

            # extra journalling post solve
            self.history.set_last_cell("mu", self.mu)
            self.history.set_last_cell("nu", nu)
            if "solution" in self.history:
                self.history.set_last_cell("solution", self.actual.step().deep_copy())
            self.history.set_last_cell("solution_norm", self.step().norm())
            self.history.set_last_cell("reduced_chi_squared", self.chi_sq())

            # test termination criteria before taking the next forward step
            if self.had_too_small_a_step():
                self.history.reason_for_termination = STEP_TOO_SMALL
                break
            if self.test_for_termination():
                self.history.reason_for_termination = TARGET_ACHIEVED
                break
            if self.test_rmsd_convergence():
                self.history.reason_for_termination = RMSD_CONVERGED
                break
            if self.n_iterations == self._max_iterations:
                self.history.reason_for_termination = MAX_ITERATIONS
                break

            h = self.step()
            expected_decrease = 0.5 * h.dot(self.mu * h - self._g)
            self.step_forward()
            self.n_iterations += 1
            self.build_up(objective_only=True)
            objective_new = self.objective()
            self.report_progress(objective_new)
            actual_decrease = self._f - objective_new
            rho = actual_decrease / expected_decrease
            if rho > 0:
                self.mu *= max(1 / 3, 1 - (2 * rho - 1) ** 3)
                nu = 2
            else:
                self.step_backward()
                self.history.del_last_row()
                if nu >= 8192:
                    self.history.reason_for_termination = MAX_TRIAL_ITERATIONS
                    break
                self.mu *= nu
                nu *= 2

            # prepare for next step
            self.build_up()

        self.calculate_esds()

    def run(self):
        self._run_core()
        self.calculate_esds()
