# Copyright (C) (2012) David Waterman, STFC Rutherford Appleton Laboratory, UK.
# This code is developed as part of the DIALS project and is provided for
# testing purposes only

from __future__ import division
from scitbx import lbfgs
from cctbx.array_family import flex
import math
import libtbx

# use lstbx classes
from scitbx.lstbx import normal_eqns, normal_eqns_solving

class refinery(object):
    '''Abstract interface for refinery objects'''

    # NOTES. A refinery is initialised with a target function. The target
    # function already contains a reflection manager (which holds the data) so
    # there's no need to pass the data in here. In fact the target function
    # class does the bulk of the work, as it also does the reflection prediction
    # to get the updated predictions on each cycle. This makes some kind of sense
    # as the target function is inextricably linked to the space in which
    # predictions are made (e.g. detector space, phi), so it is not general
    # enough to sit abstractly above the prediction.

    # This keeps the refinery simple and able to be focused only on generic
    # features of managing a refinement run, like reporting results and checking
    # termination criteria

    # The prediction values come from a prediction_parameterisation object.
    # This is also referred to by the target function, but it makes sense for
    # refinery to be able to refer to it directly. So refinery should keep a
    # separate link to its prediction_parameterisation

    def __init__(self, target, prediction_parameterisation, log=None,
                 verbosity = 0):

        self._parameters = prediction_parameterisation
        self._target = target
        self.x = flex.double(self._parameters.get_p())
        self._log = log
        self._step = 0
        self._verbosity = verbosity
        self._target_achieved = False
        self.compute_functional_and_gradients()
        if self._verbosity > 0.: self.print_table_row()

    def get_num_steps(self):
        return self._step

    def compute_functional_and_gradients(self):

        # set current parameter values
        self._parameters.set_p(self.x)

        # do reflection prediction
        self._target.predict()

        # compute target function and its gradients
        self._f, self._g = self._target.compute_functional_and_gradients()

        if self._verbosity > 1:
            print "Function evaluation"
            msg = "  Params: " + "%.5f " * len(self._parameters)
            print msg % tuple(self._parameters.get_p())
            print "  L: %.5f" % self._f
            msg = "  dL/dp: " + "%.5f " * len(tuple(self._g))
            print msg % tuple(self._g)

        return self._f, flex.double(self._g)

    def callback_after_step(self, minimizer):
        '''returns True to terminate the refinement'''

        self._step += 1
        if self._verbosity > 0.: self.print_table_row()

        # delegate this to the target class
        self._target_achieved = self._target.achieved()
        return self._target_achieved

    def print_table_row(self):
        '''print useful output in the form of a tab separated table'''

        params = self._parameters.get_p()
        if self._step == 0: # first call, print the column headings
            header = ("Step Nref Residual RMSD_X RMSD_Y RMSD_phi " +
                     "Param_%02d " * len(params))
            print header % tuple(range(1, len(params) + 1))

        line = ("%d " + "%d " + "%.5f " + "%.5f " * 3 +
                "%.5f " * len(params))

        rmsds = self._target.rmsds()
        print line % tuple((self._step,
                           self._target.get_num_reflections(),
                           self._f,
                           rmsds[0],
                           rmsds[1],
                           rmsds[2]) + tuple(params))

    def run(self):
        '''To be implemented by derived class'''

        # Specify a minimizer and its parameters, and run
        raise RuntimeError, "implement me"

class simple_lbfgs(refinery):
    '''Refinery implementation, using cctbx LBFGS with basic settings'''

    def run(self):

        #TODO convert file file handling lines to use of 'with'?
        ref_log = None
        if self._log: ref_log = open(self._log, "w")
        self.minimizer = lbfgs.run(target_evaluator=self, log=ref_log)
        if self._log: ref_log.close()

class lbfgs_curvs(refinery):
    '''LBFGS refinery using curvatures'''

    def run(self):

        #TODO convert file file handling lines to use of 'with'?
        ref_log = None
        if self._log: ref_log = open(self._log, "w")
        self.diag_mode = "always"
        self.minimizer = lbfgs.run(target_evaluator=self, log=ref_log)
        if self._log: ref_log.close()

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


class adapt_lstbx(
    refinery,
    normal_eqns.non_linear_ls,
    normal_eqns.non_linear_ls_mixin):
    '''refinery using lstbx
    testing the use of lstbx methods for solving the normal equations

    Here following the example of lstbx.tests.test_problems.exponential_fit
    to interface to lstbx code.

    There is much duplication of code with other classes while I figure
    out what the interface and inheritance hierarchy should be'''

    def __init__(self, target, prediction_parameterisation, log=None,
                 verbosity = 0):

        refinery.__init__(self, target, prediction_parameterisation, log=None,
                 verbosity = 0)

        # required for restart to work (do I need that method?)
        self.x_0 = self.x.deep_copy()

        normal_eqns.non_linear_ls.__init__(self, n_parameters = len(self._parameters))

        # determine overall scale factor so that the objective is approx in
        # [0,1]
        self._target.predict()
        self._scale = 1./math.sqrt(self._target.compute_functional_and_gradients()[0])


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
        self._parameters.set_p(self.x)

        # do reflection prediction
        self._target.predict()

        if self._verbosity > 1:
            print "Function evaluation"
            msg = "  Params: " + "%.5f " * len(self._parameters)
            print msg % tuple(self._parameters.get_p())
            print "  L: %.5f" % self._f
            msg = "  dL/dp: " + "%.5f " * len(tuple(self._g))
            print msg % tuple(self._g)

        residuals, jacobian, weights = \
            self._target.compute_residuals_and_gradients()

        #print "sum of residuals", sum(residuals)
        #print "objective", 0.5* sum(weights * residuals**2)
        #print "unweighted objective", 0.5* sum(residuals**2)
        #print "scaled objective", 0.5* sum(weights * self._scale**2 * residuals**2)
        #print "scale factor", self._scale

        # apply overall scale factor
        residuals *= self._scale
        jacobian *= self._scale
        #weights *= (self._scale)**2.

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

class gn_iterations(adapt_lstbx, normal_eqns_solving.iterations):

    track_step = False
    track_gradient = False
    track_all = False
    n_max_iterations = 100
    gradient_threshold = 1.e-10
    step_threshold = None
    damping_value = 0.0007
    max_shift_over_esd = 15
    convergence_as_shift_over_esd = 1e-5

    # override the base class __init__ as I don't want to start
    # refinement on initialisation
    def __init__(self, target, prediction_parameterisation, log=None,
                 verbosity = 0, **kwds):
        """
        """

        adapt_lstbx.__init__(self, target, prediction_parameterisation, log=None,
                 verbosity = 0)

        libtbx.adopt_optional_init_args(self, kwds)
        if self.track_all: self.track_step = self.track_gradient = True
        self.non_linear_ls = normal_eqns_solving.journaled_non_linear_ls(
            non_linear_ls = self, journal = self, track_gradient = self.track_gradient,
            track_step = self.track_step)

    def run(self):
        self.n_iterations = 0
        while self.n_iterations < self.n_max_iterations:
            self.non_linear_ls.build_up()
            if self.has_gradient_converged_to_zero():
                print "Gradient converged to zero"
                break
            if self.non_linear_ls.callback_after_step(None):
                print "RMSD target achieved"
                break
            self.non_linear_ls.solve()
            if self.had_too_small_a_step(): break
            self.non_linear_ls.step_forward()
            self.n_iterations += 1

    def __str__(self):
        return "pure Gauss-Newton"
