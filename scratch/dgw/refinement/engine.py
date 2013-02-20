# Copyright (C) (2012) David Waterman, STFC Rutherford Appleton Laboratory, UK.
# This code is developed as part of the DIALS project and is provided for
# testing purposes only

from __future__ import division
from scitbx import lbfgs
from cctbx.array_family import flex

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

        self.minimizer = lbfgs.run(target_evaluator=self, log=self._log)

class lbfgs_curvs(refinery):
    '''LBFGS refinery using curvatures'''

    def run(self):

        self.diag_mode = "always"
        self.minimizer = lbfgs.run(target_evaluator=self, log=self._log)

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
