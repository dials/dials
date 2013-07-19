#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

from __future__ import division
from scitbx import lbfgs
from cctbx.array_family import flex
import libtbx

# use lstbx classes
from scitbx.lstbx import normal_eqns, normal_eqns_solving

class Refinery(object):
    '''Abstract interface for Refinery objects'''

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

    def __init__(self, target, prediction_parameterisation, log = None,
                 verbosity = 0, track_step = False,
                 track_gradient = False):

        # reference to PredictionParameterisation and Target objects
        self._parameters = prediction_parameterisation
        self._target = target

        # initial parameter values
        self.x = flex.double(self._parameters.get_p())

        # undefined initial functional and gradients values
        self._f = None
        self._g = None

        # filename for an optional log file
        self._log = log

        self._verbosity = verbosity

        self._target_achieved = False

        # attributes for journalling functionality, based on lstbx's
        # journaled_non_linear_ls class
        self._step = 0
        self.num_reflections_history = []
        self.objective_history = flex.double()
        self.gradient_history = [] if track_gradient else None
        self.gradient_norm_history = flex.double()
        self.step_history = [] if track_step else None
        self.step_norm_history = flex.double()
        self.parameter_vector_history = []
        self.parameter_vector_norm_history = flex.double()
        self.rmsd_history = []

        self.prepare_for_step()

    def get_num_steps(self):
        return self._step

    def prepare_for_step(self):
        '''Update the parameterisation and prepare the target function'''

        # set current parameter values
        self._parameters.set_p(self.x)

        # do reflection prediction
        self._target.predict()

    def update_journal(self):
        '''Append latest step information to the journal attributes'''

        # add step quantities to journal
        self._step += 1
        self.num_reflections_history.append(self._target.get_num_reflections())
        self.rmsd_history.append(self._target.rmsds())
        self.parameter_vector_history.append(self._parameters.get_p())
        self.objective_history.append(self._f)
        if self.gradient_history is not None:
            self.gradient_history.append(self._g)

    def test_for_termination(self):
        '''Return True if refinement should be terminated'''

        # Basic version delegate to the Target class. Derived classes may
        # implement other termination criteria
        self._target_achieved = self._target.achieved()

        return self._target_achieved

    def test_rmsd_convergence(self):
        '''Test for convergence of RMSDs'''

        # http://en.wikipedia.org/wiki/
        # Non-linear_least_squares#Convergence_criteria
        try:
            r1 = self.rmsd_history[-1]
            r2 = self.rmsd_history[-2]
        except IndexError:
            return False

        tests = [abs((e[1] - e[0])/e[1])  < 0.0001 for e in zip(r1, r2)]

        return all(tests)

    def print_step(self):
        '''print information about the current step'''

        print "Function evaluation"
        msg = "  Params: " + "%.5f " * len(self._parameters)
        print msg % tuple(self._parameters.get_p())
        print "  L: %.5f" % self._f
        msg = "  dL/dp: " + "%.5f " * len(tuple(self._g))
        print msg % tuple(self._g)

    def print_table(self):
        '''print useful output in the form of a space-separated table'''

        print
        print "Refinement steps"
        print "----------------"
        header = ("Step Nref Residual RMSD_X RMSD_Y RMSD_phi " +
                     "Param_%02d " * len(self._parameters))
        print header % tuple(range(1, len(self._parameters) + 1))

        for i in range(self._step):
            dat = (i + 1,) + (self.num_reflections_history[i],) + \
                  (self.objective_history[i],) + \
                  tuple(self.rmsd_history[i]) + \
                  tuple(self.parameter_vector_history[i])
            print  ("%d " + "%d " + "%.5f " + "%.5f " * 3 +
                "%.5f " * len(self._parameters)) % dat

    def run(self):
        '''
        To be implemented by derived class. It is expected that each step of
        refinement be preceeded by a call to prepare_for_step and followed by
        calls to update_journal and test_for_termination (in that order).
        '''

        # Specify a minimizer and its parameters, and run
        raise RuntimeError, "implement me"

class AdaptLbfgs(Refinery):
    '''Adapt Refinery for L-BFGS minimiser'''

    def compute_functional_and_gradients(self):

        self.prepare_for_step()

        # compute target function and its gradients
        self._f, self._g = self._target.compute_functional_and_gradients()

        if self._verbosity > 1: self.print_step()

        return self._f, flex.double(self._g)

    def callback_after_step(self, minimizer):
        '''
        Do journalling, evaluate rmsds and return True if the target is
        reached to terminate the refinement.
        '''

        self.update_journal()
        if self._verbosity > 0.:
            print "Step", self._step

        return self.test_for_termination()

class SimpleLBFGS(AdaptLbfgs):
    '''Refinery implementation, using cctbx LBFGS with basic settings'''

    def run(self):

        #TODO convert file file handling lines to use of 'with'?
        ref_log = None
        if self._log: ref_log = open(self._log, "w")
        self.minimizer = lbfgs.run(target_evaluator=self, log=ref_log)
        if self._log: ref_log.close()

        if self._verbosity > 0.: self.print_table()

        return

class LBFGScurvs(AdaptLbfgs):
    '''Refinery implementation using cctbx LBFGS with curvatures'''

    def run(self):

        #TODO convert file file handling lines to use of 'with'?
        ref_log = None
        if self._log: ref_log = open(self._log, "w")
        self.diag_mode = "always"
        self.minimizer = lbfgs.run(target_evaluator=self, log=ref_log)
        if self._log: ref_log.close()

        if self._verbosity > 0.: self.print_table()

        return

    def compute_functional_gradients_diag(self):

        f, g = self.compute_functional_and_gradients()
        curvs = self.curvatures()

        diags = 1. / curvs

        if self._verbosity > 1:
            msg = "  curv: " +  "%.5f " * len(tuple(curvs))
            print msg % tuple(curvs)

        return self._f, flex.double(self._g), diags

    def curvatures(self):

        # This relies on compute_functional_and_gradients being called first
        # (in order to set the parameters and update predictions).
        return(flex.double(self._target.curvatures()))


class AdaptLstbx(
    Refinery,
    normal_eqns.non_linear_ls,
    normal_eqns.non_linear_ls_mixin):
    '''Adapt Refinery for lstbx'''

    def __init__(self, target, prediction_parameterisation, log=None,
                 verbosity = 0, track_step = False,
                 track_gradient = False):

        Refinery.__init__(self, target, prediction_parameterisation,
                 log=None, verbosity = verbosity, track_step = False,
                 track_gradient = False)

        # required for restart to work (do I need that method?)
        self.x_0 = self.x.deep_copy()

        normal_eqns.non_linear_ls.__init__(self, n_parameters = len(self._parameters))

        # determine overall scale factor for the gradient threshold
        self._target.predict()
        self._scale = self._target.compute_functional_and_gradients()[0]


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

        # get calculations from the target
        residuals, jacobian, weights = \
            self._target.compute_residuals_and_gradients()

        #Reset the state to construction time, i.e. no equations accumulated
        self.reset()

        if objective_only:
            self.add_residuals(residuals, weights)
        else:
            self.add_equations(residuals, jacobian, weights)

    def step_forward(self):
        self.old_x = self.x.deep_copy()
        self.x += self.step()

    def step_backward(self):
        assert self.old_x is not None
        self.x, self.old_x = self.old_x, None
        self._step -= 1

class GaussNewtonIterations(AdaptLstbx, normal_eqns_solving.iterations):
    '''Refinery implementation, using lstbx Gauss Newton iterations'''

    n_max_iterations = 100
    gradient_threshold = 1.e-10
    step_threshold = None
    damping_value = 0.0007
    max_shift_over_esd = 15
    convergence_as_shift_over_esd = 1e-5
    reduced_chi_squared_history = flex.double()

    def __init__(self, target, prediction_parameterisation, log=None,
                 verbosity = 0, track_step = False,
                 track_gradient = False, **kwds):

        AdaptLstbx.__init__(self, target, prediction_parameterisation,
                 log = None, verbosity = verbosity, track_step = False,
                 track_gradient = False)

        # scale gradient threshold for this problem
        self.gradient_threshold *= self._scale

        libtbx.adopt_optional_init_args(self, kwds)

    def run(self):
        self.n_iterations = 0
        while self.n_iterations < self.n_max_iterations:
            self.build_up()

            # set functional and gradients for the step (to be added
            # to the history at callback_after_step())
            self._f = self.objective()
            self._g = -self.opposite_of_gradient()

            # journalling prior to solve
            self.parameter_vector_norm_history.append(
              self.parameter_vector_norm())
            self.gradient_norm_history.append(
              self.opposite_of_gradient().norm_inf())

            if self._verbosity > 2:
                print "The normal matrix for the current step is:"
                print self.normal_matrix_packed_u().\
                    matrix_packed_u_as_symmetric().\
                    as_scitbx_matrix().matlab_form(format=None,
                    one_row_per_line=True)
                print

            self.update_journal()
            if self._verbosity > 0.: print "Step", self._step
            if self._verbosity > 1: self.print_step()

            # solve the normal equations
            self.solve()

            # extra journalling post solve
            if self.step_history is not None:
              self.step_history.append(self.actual.step().deep_copy())
            self.step_norm_history.append(self.step().norm())
            self.reduced_chi_squared_history.append(self.chi_sq())

            # test termination criteria
            if self.test_for_termination():
                reason_for_termination = "RMSD target achieved"
                break
            if self.test_rmsd_convergence():
                reason_for_termination = "RMSD no longer decreasing"
                break
            #if self.has_gradient_converged_to_zero():
            #    print "Gradient converged to zero"
            #    break
            if self.had_too_small_a_step():
                reason_for_termination = "Step too small"
                break

            # prepare for next step
            self.step_forward()
            self.n_iterations += 1

        # print output table
        if self._verbosity > 0.:
            self.print_table()
            print reason_for_termination

        # invert normal matrix from N^-1 = (U^-1)(U^-1)^T
        cf = self.step_equations().cholesky_factor_packed_u()
        cf_inv = cf.matrix_packed_u_as_upper_triangle().\
            matrix_inversion()
        nm_inv = cf_inv.matrix_multiply_transpose(cf_inv)

        # keep the estimated parameter variance-covariance matrix
        self.parameter_var_cov = \
            self.reduced_chi_squared_history[-1] * nm_inv

        return
