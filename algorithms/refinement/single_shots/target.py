#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

# python and cctbx imports
from __future__ import division
from scitbx import matrix
from math import pi, sqrt
from cctbx.array_family import flex
import random

# dials imports
from dials.algorithms.spot_prediction import ray_intersection
from dials.algorithms.refinement.target import Target
from dials.algorithms.refinement.target import ReflectionManager
from dials.algorithms.refinement.target import ObsPredMatch
from dials.algorithms.refinement.target import ObservationPrediction

# constants
TWO_PI = 2.0 * pi

class LeastSquaresXYResidualWithRmsdCutoff(Target):
    """An implementation of the target class providing a least squares residual
    in terms of detector impact position X, Y only, terminating on achieved
    rmsd (or on intrisic convergence of the chosen minimiser)"""

    rmsd_names = ["RMSD_X", "RMSD_Y"]

    def __init__(self, reflection_predictor, detector, ref_man,
                 prediction_parameterisation,
                 frac_binsize_cutoff=0.33333,
                 absolute_cutoffs=None):

        Target.__init__(self, reflection_predictor, detector, ref_man,
                        prediction_parameterisation)

        # Set up the RMSD achieved criterion
        if not absolute_cutoffs:
            pixel_sizes = [p.get_pixel_size() for p in detector]
            min_px_size_x = min(e[0] for e in pixel_sizes)
            min_px_size_y = min(e[1] for e in pixel_sizes)
            self._binsize_cutoffs = [min_px_size_x * frac_binsize_cutoff,
                                     min_px_size_y * frac_binsize_cutoff]
        else:
            self._binsize_cutoffs = absolute_cutoffs[:2]

        # Quantities to cache each step
        self._rmsds = None
        self._matches = None

    def compute_residuals_and_gradients(self):
        """return the vector of residuals plus their gradients
        and weights for non-linear least squares methods"""

        self._matches = self._reflection_manager.get_matches()

        # return residuals and weights as 1d flex.double vectors
        nelem = len(self._matches) * 2
        residuals = flex.double(nelem)
        jacobian_t = flex.double(flex.grid(
            len(self._prediction_parameterisation), nelem))
        weights = flex.double(nelem)

        for i, m in enumerate(self._matches):
            residuals[2*i] = m.Xresid
            residuals[2*i + 1] = m.Yresid
            #residuals[3*i + 2] = m.Phiresid

            # are these the right weights? Or inverse, or sqrt?
            weights[2*i] = m.weightXo
            weights[2*i + 1] = m.weightYo
            #weights[3*i + 2] = m.weightPhio

            # m.gradients is a nparam length list, each element of which is a
            # doublet of values, (dX/dp_n, dY/dp_n)
            dX_dp, dY_dp = zip(*m.gradients)

            # FIXME Here we paste columns into the Jacobian transpose then take
            # its transpose when done. This seems inefficient: can we just start
            # with the Jacobian and fill elements sequentially, using row-major
            # order to ensure the values are filled in the right order?

            # fill jacobian elements here.
            jacobian_t.matrix_paste_column_in_place(flex.double(dX_dp), 2*i)
            jacobian_t.matrix_paste_column_in_place(flex.double(dY_dp), 2*i+1)
            #jacobian_t.matrix_paste_column_in_place(flex.double(dPhi_dp), 3*i+2)

        # We return the Jacobian, not its transpose.
        jacobian_t.matrix_transpose_in_place()

        return(residuals, jacobian_t, weights)

    def compute_functional_and_gradients(self):
        """calculate the value of the target function and its gradients"""

        self._matches = self._reflection_manager.get_matches()
        self._nref = self.get_num_reflections()

        # This is a hack for the case where nref=0. This should not be necessary
        # if bounds are provided for parameters to stop the algorithm exploring
        # unreasonable regions of parameter space where no predictions exist.
        # Unfortunately the L-BFGS line search does make such extreme trials.
        if self._nref == 0:
            return 1.e12, [1.] * len(self._prediction_parameterisation)

        # compute target function
        L = 0.5 * sum([m.weightXo * m.Xresid2 +
                       m.weightYo * m.Yresid2
                       for m in self._matches])

        # prepare list of gradients
        dL_dp = [0.] * len(self._prediction_parameterisation)

        # the gradients wrt each parameter are stored with the matches
        for m in self._matches:

            for j, (grad_X, grad_Y) in enumerate(m.gradients):
                dL_dp[j] += (m.weightXo * m.Xresid * grad_X +
                             m.weightYo * m.Yresid * grad_Y)

        return (L, dL_dp)

    def curvatures(self):
        """First order approximation to the diagonal of the Hessian based on the
        least squares form of the target"""

        # This is a hack for the case where nref=0. This should not be necessary
        # if bounds are provided for parameters to stop the algorithm exploring
        # unreasonable regions of parameter space where there are no predictions
        if self._nref == 0:
            return [1.] * len(self._prediction_parameterisation)

        # prepare lists of gradients and curvatures
        curv = [0.] * len(self._prediction_parameterisation)

        # for each reflection, get the approximate curvatures wrt each parameter
        for m in self._matches:

            for j, (grad_X, grad_Y) in enumerate(m.gradients):
                curv[j] += (m.weightXo * grad_X**2 +
                            m.weightYo * grad_Y**2)

        # Curvatures of zero will cause a crash, because their inverse is taken.
        assert all([c > 0.0 for c in curv])

        return curv

    def rmsds(self):
        """calculate unweighted RMSDs"""

        if not self._matches:
            self._matches = self._reflection_manager.get_matches()

        n = self._reflection_manager.get_accepted_reflection_count()

        resid_x = sum((m.Xresid2 for m in self._matches))
        resid_y = sum((m.Yresid2 for m in self._matches))

        # cache rmsd calculation for achieved test
        self._rmsds = (sqrt(resid_x / n),
                       sqrt(resid_y / n))

        return self._rmsds

    def achieved(self):
        """RMSD criterion for target achieved """
        r = self._rmsds if self._rmsds else self.rmsds()

        # reset cached rmsds to avoid getting out of step
        self._rmsds = None

        if (r[0] < self._binsize_cutoffs[0] and
            r[1] < self._binsize_cutoffs[1]):
            return True
        return False

class ReflectionManagerXY(ReflectionManager):
    """Overloads for a Reflection Manager that does not exclude
    reflections too close to the spindle, and reports only information
    about X, Y residuals"""

    def _spindle_beam_plane_normal(self):
        """There is no goniometer, so overload to return None"""

        return None

    def _id_refs_to_keep(self, obs_data):
        """For this version of the class, do nothing. We don't want to
        exclude reflections close to the spindle, as the spindle may
        not exist"""

        inc = [i for i, ref in enumerate(obs_data)]

        return inc

    def _create_working_set(self, indices):
        """Make a subset of the indices of reflections to use in refinement.

        This version ignores nref_per_degree"""

        working_indices = indices
        sample_size = len(working_indices)

        # set maximum sample size
        if self._max_num_obs:
            if sample_size > self._max_num_obs:
                sample_size = self._max_num_obs

        # sample the data and record the sample size
        if sample_size < len(working_indices):
            self._sample_size = sample_size
            working_indices = random.sample(working_indices,
                                            self._sample_size)
        return(working_indices)

    def get_matches(self, silent = False):
        """For every observation matched with a prediction return all data"""

        l = [obs for v in self._obs_pred_pairs.values() for obs in v.obs if obs.is_matched]

        if self._verbosity > 2 and len(l) > 20 and not silent:
            print "Listing predictions matched with observations for " + \
                  "the first 20 reflections:"
            print "H, K, L, Xresid, Yresid, weightXo, weightYo"
            fmt = "(%3d, %3d, %3d) %5.3f %5.3f %5.3f %5.3f"
            for i in xrange(20):
                e = l[i]
                msg = fmt % tuple(e.H + (e.Xresid,
                                 e.Yresid,
                                 e.weightXo,
                                 e.weightYo))
                print msg
            print

        return l
