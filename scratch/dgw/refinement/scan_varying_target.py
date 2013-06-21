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

# constants
TWO_PI = 2.0 * pi

class Target(object):
    '''Abstract interface for a target function class

    A Target object will be used by a Refinery. It will refer to a Reflection
    Manager to get a list of observations. It will perform reflection prediction
    on those observations and update the reflection manager with those
    predictions. It will then query the reflection manager to get a list of
    accepted reflections with observed and calculated positions. These are the
    reflections for use in refinement. It obtains the gradients of reflection
    positions from a relevant prediction parameterisation object. With all of
    this information in place, it calculates the value of the target function,
    the gradients of the target function and auxiliary information (e.g. RMSDs).

    Concrete instances of this class implement the actual target function
    calculation. The base class should not determine what type of target
    function is used (e.g. least squares target), or limit whether the object
    is used for a detector space / phi residual, or a reciprocal space residual.
    This should all be set by a derived class.
    '''

    def __init__(self, reflection_predictor, detector,
                 ref_manager,
                 prediction_parameterisation):

        self._reflection_predictor = reflection_predictor
        self._detector = detector
        self._H = ref_manager
        self._prediction_parameterisation = prediction_parameterisation

    def predict(self):
        '''perform reflection prediction and update the reflection manager'''

        # update the reflection_predictor and the prediction parameterisation
        # with the scan-independent part of the current geometry
        self._reflection_predictor.update()
        self._prediction_parameterisation.prepare()

        # reset the 'use' flag for all observations
        self._H.reset_accepted_reflections()

        # loop over all reflections in the manager
        for h in self._H.get_indices():

            # loop over observations of this hkl
            for obs in self._H.get_obs(h):

                # get the observation image number
                frame = obs.frame_o

                # duck-typing for scan varying version of
                # prediction_parameterisation
                try:

                    # compose the prediction parameterisation at the requested
                    # image number
                    self._prediction_parameterisation.compose(frame)

                    # extract UB matrix
                    UB = self._prediction_parameterisation.get_UB(frame)

                    # predict for this hkl at setting UB
                    predictions = self._reflection_predictor.predict(h, UB)

                except AttributeError:

                    # predict for this hkl
                    predictions = self._reflection_predictor.predict(h)

                # obtain the impact positions, currently assuming reflections
                # only intersect panel 0
                impacts = ray_intersection(self._detector, predictions, panel=0)

                # find the prediction with the right 'entering' flag
                try:
                    i = [x.entering == obs.entering for x in impacts].index(True)
                except ValueError:
                    # we don't have a prediction for this obs
                    continue

                ref = impacts[i]
                Xc, Yc = ref.image_coord_mm

                # do not wrap around multiples of 2*pi; keep the full rotation
                # from zero to differentiate repeat observations.
                resid = ref.rotation_angle - (obs.Phio % TWO_PI)
                Phic = obs.Phio + resid
                Sc = matrix.col(ref.beam_vector)

                # calculate gradients for this reflection
                grads = self._prediction_parameterisation.get_gradients(
                                                h, Sc, Phic, frame)

                # store all this information in the matched obs-pred pair
                obs.update_prediction(Xc, Yc, Phic, Sc, grads)

        if self._H.first_update:

            # delete all obs-pred pairs from the manager that do not
            # have a prediction
            self._H.strip_unmatched_observations()

            self._H.first_update = False

        return

    def get_num_reflections(self):
        '''return the number of reflections currently used in the calculation'''

        # delegate to the reflection manager
        return self._H.get_accepted_reflection_count()

    def compute_functional_and_gradients(self):
        '''calculate the target function value and its gradients'''

        # To be implemented by a derived class
        raise RuntimeError, 'implement me'

    def achieved(self):
        '''return True to terminate the refinement. To be implemented by
        a derived class'''

        return False

class LeastSquaresPositionalResidualWithRmsdCutoff(Target):
    '''An implementation of the target class providing a least squares residual
    in terms of detector impact position X, Y and phi, terminating on achieved
    rmsd (or on intrisic convergence of the chosen minimiser)'''

    def __init__(self, reflection_predictor, detector, ref_man,
                 prediction_parameterisation,
                 image_width):

        Target.__init__(self, reflection_predictor, detector, ref_man,
                        prediction_parameterisation)

        # Scale units for the RMSD achieved criterion
        self._pixel_size = detector.get_pixel_size()
        self._image_width = image_width

    def compute_residuals_and_gradients(self):
        '''return the vector of residuals plus their gradients
        and weights for non-linear least squares methods'''

        self._matches = self._H.get_matches()

        # return residuals and weights as 1d flex.double vectors
        nelem = len(self._matches) * 3
        residuals = flex.double(nelem)
        jacobian_t = flex.double(flex.grid(
            len(self._prediction_parameterisation), nelem))
        weights = flex.double(nelem)

        for i, m in enumerate(self._matches):
            residuals[3*i] = m.Xresid
            residuals[3*i + 1] = m.Yresid
            residuals[3*i + 2] = m.Phiresid

            # are these the right weights? Or inverse, or sqrt?
            weights[3*i] = m.weightXo
            weights[3*i + 1] = m.weightYo
            weights[3*i + 2] = m.weightPhio

            # m.gradients is a nparam length list, each element of which is a
            # triplet of values, (dX/dp_n, dY/dp_n, dPhi/dp_n)
            dX_dp, dY_dp, dPhi_dp = zip(*m.gradients)

            # FIXME Here we paste columns into the Jacobian transpose then take
            # its transpose when done. This seems inefficient: can we just start
            # with the Jacobian and fill elements sequentially, using row-major
            # order to ensure the values are filled in the right order?

            # fill jacobian elements here.
            jacobian_t.matrix_paste_column_in_place(flex.double(dX_dp), 3*i)
            jacobian_t.matrix_paste_column_in_place(flex.double(dY_dp), 3*i + 1)
            jacobian_t.matrix_paste_column_in_place(flex.double(dPhi_dp), 3*i + 2)

        # We return the Jacobian, not its transpose.
        jacobian_t.matrix_transpose_in_place()

        return(residuals, jacobian_t, weights)

    def compute_functional_and_gradients(self):
        '''calculate the value of the target function and its gradients'''

        self._matches = self._H.get_matches()
        self._nref = self.get_num_reflections()

        # This is a hack for the case where nref=0. This should not be necessary
        # if bounds are provided for parameters to stop the algorithm exploring
        # unreasonable regions of parameter space where no predictions exist.
        # Unfortunately the L-BFGS line search does make such extreme trials.
        if self._nref == 0:
            return 1.e12, [1.] * len(self._prediction_parameterisation)

        # compute target function
        L = 0.5 * sum([m.weightXo * m.Xresid2 +
                       m.weightYo * m.Yresid2 +
                       m.weightPhio * m.Phiresid2
                       for m in self._matches])

        # prepare list of gradients
        dL_dp = [0.] * len(self._prediction_parameterisation)

        # the gradients wrt each parameter are stored with the matches
        for m in self._matches:

            for j, (grad_X, grad_Y, grad_Phi) in enumerate(m.gradients):
                dL_dp[j] += (m.weightXo * m.Xresid * grad_X +
                             m.weightYo * m.Yresid * grad_Y +
                             m.weightPhio * m.Phiresid * grad_Phi)

        return (L, dL_dp)

    def curvatures(self):
        '''First order approximation to the diagonal of the Hessian based on the
        least squares form of the target'''

        # This is a hack for the case where nref=0. This should not be necessary
        # if bounds are provided for parameters to stop the algorithm exploring
        # unreasonable regions of parameter space where there are no predictions
        if self._nref == 0:
            return [1.] * len(self._prediction_parameterisation)

        # prepare lists of gradients and curvatures
        curv = [0.] * len(self._prediction_parameterisation)

        # for each reflection, get the approximate curvatures wrt each parameter
        for m in self._matches:

            for j, (grad_X, grad_Y, grad_Phi) in enumerate(m.gradients):
                curv[j] += (m.weightXo * grad_X**2 +
                            m.weightYo * grad_Y**2 +
                            m.weightPhio * grad_Phi**2)

        # Curvatures of zero will cause a crash, because their inverse is taken.
        assert all([c > 0.0 for c in curv])

        return curv

    def rmsds(self):
        '''calculate unweighted RMSDs'''

        n = self._H.get_accepted_reflection_count()

        resid_x = sum((m.Xresid2 for m in self._matches))
        resid_y = sum((m.Yresid2 for m in self._matches))
        resid_phi = sum((m.Phiresid2 for m in self._matches))

        return sqrt(resid_x / n), sqrt(resid_y / n), sqrt(resid_phi / n)

    def achieved(self):
        '''RMSD criterion for target achieved '''
        r = self.rmsds()
        if (r[0] < self._pixel_size[0] * 0.33333 and
            r[1] < self._pixel_size[1] * 0.33333 and
            r[2] < self._image_width * 0.33333):
            return True
        return False

class ObsPredMatch:
    '''
    A bucket class containing data for a prediction that has been
    matched to an observation.

    This contains all the raw material needed to calculate the target function
    value and gradients
    '''

    # initialise with an observation
    def __init__(self, hkl, entering, frame, Xo, sigXo, weightXo,
                            Yo, sigYo, weightYo,
                            Phio, sigPhio, weightPhio):
        self.H = hkl
        self.entering = entering
        self.frame_o = frame
        self.Xo = Xo
        self.sigXo = sigXo
        self.weightXo = weightXo

        self.Yo = Yo
        self.sigYo = sigYo
        self.weightYo = weightYo

        self.Phio = Phio
        self.sigPhio = sigPhio
        self.weightPhio = weightPhio

        self.Xc = None
        self.Yc = None
        self.Phic = None
        self.Sc = None

        # gradients will be a list, of length equal to the number of free
        # parameters, whose elements are triplets (dX/dp, dY/dp, dPhi/dp)
        self.gradients = None

        self.Yresid = None
        self.Yresid2 = None
        self.Xresid = None
        self.Xresid2 = None
        self.Phiresid = None
        self.Phiresid2 = None

        self.is_matched = False

    # update with a prediction
    def update_prediction(self, Xc, Yc, Phic, Sc, gradients):

        self.Xc = Xc
        self.Yc = Yc
        self.Phic = Phic
        self.Sc = Sc

        self.gradients = gradients

        # calculate residuals
        self.Xresid = Xc - self.Xo
        self.Xresid2 = self.Xresid**2
        self.Yresid = Yc - self.Yo
        self.Yresid2 = self.Yresid**2
        self.Phiresid = Phic - self.Phio
        self.Phiresid2 = self.Phiresid**2

        self.is_matched = True

    def reset(self):

        '''Flag this observation to not be used'''
        self.is_matched = False

class ObservationPrediction(object):
    '''A helper class for the reflection manager to contain information about
    the unique observations of particular hkl paired with the currently
    predicted values.

    Reflections are either exiting or entering the Ewald sphere. This is
    calculated by the reflection manager and passed in, and used here to
    identify which observations to update with a new prediction'''

    def __init__(self, H, entering, frame,
                       Xo, sigXo, weightXo,
                       Yo, sigYo, weightYo,
                       Phio, sigPhio, weightPhio):
        '''initialise from one observation'''
        assert isinstance(entering, bool)
        self.H = H

        self.obs = [ObsPredMatch(H, entering, frame,
                                 Xo, sigXo, weightXo,
                                 Yo, sigYo, weightYo,
                                 Phio, sigPhio, weightPhio)]

    def __len__(self):

        return len(self.obs)

    def __iter__(self):

        # make this class iterable through the observations
        return iter(self.obs)

    def get_num_pairs(self):
        '''Count the number of observations that are paired with a prediction'''
        return sum(1 for x in self.obs if x.is_matched)

    def reset_predictions(self):
        '''Set the use flag to false for all observations, so that after
        updating, any observations that still do not have a prediction are
        flagged to be removed from calculation of residual and gradients.'''

        map(lambda x: x.reset(), self.obs)

    def remove_unmatched_obs(self):
        '''Remove observations without a matching prediction'''

        self.obs = [x for x in self.obs if x.is_matched]

    def add_observation(self, entering, frame, Xo, sigXo, weightXo,
                              Yo, sigYo, weightYo,
                              Phio, sigPhio, weightPhio):
        '''Add another observation of this reflection'''

        self.obs.append(ObsPredMatch(self.H, entering, frame,
                                     Xo, sigXo, weightXo,
                                     Yo, sigYo, weightYo,
                                     Phio, sigPhio, weightPhio))


class ReflectionManager(object):
    '''A class to maintain information about observed and predicted
    reflections for refinement.'''

    def __init__(self, Ho, entering_o, Frameo, So,
                       Xo, sigXo,
                       Yo, sigYo,
                       Phio, sigPhio,
                       beam, gonio, scan,
                       verbosity=0, nref_per_degree = None):

        # check the observed values
        Ho = list(Ho)
        So = list(So)
        Xo = list(Xo)
        sigXo = list(sigXo)
        Yo = list(Yo)
        sigYo = list(sigYo)
        Phio = list(Phio)
        sigPhio = list(sigPhio)
        Frameo = list(Frameo)
        entering_o = list(entering_o)
        assert(len(So) == \
               len(Xo) == \
               len(sigXo) == \
               len(Yo) == \
               len(sigYo) == \
               len(Phio) == \
               len(sigPhio) == \
               len(Frameo) == \
               len(Ho))

        # track whether this is the first update of predictions or not
        self.first_update = True

        # set verbosity
        self._verbosity = verbosity

        # keep references to the beam, goniometer and scan models (for
        # reflection exclusion and subsetting)
        self._beam = beam
        self._gonio = gonio
        self._scan = scan

        # find vector normal to the spindle-beam plane for the initial model
        self._vecn = self._spindle_beam_plane_normal()

        # exclude reflections that fail inclusion criteria
        obs_data = zip(Ho, So, Xo, sigXo, Yo, sigYo, Phio, sigPhio)
        self._obs_data = self._remove_excluded_obs(obs_data)
        self._sample_size = len(self._obs_data)

        # choose a random subset of data for refinement
        (Ho, So, Xo, sigXo, Yo, sigYo, Phio, sigPhio) = \
                zip(*self._create_working_set(nref_per_degree))

        # store observation information in a dict of observation-prediction
        # pairs (prediction information will go in here later)
        self._H = {}
        for i, h in enumerate(Ho):
            entering = So[i].dot(self._vecn) < 0.
            if h not in self._H:
                self._H[h] = ObservationPrediction(h, entering, Frameo[i],
                                        Xo[i], sigXo[i], 1./sigXo[i]**2,
                                        Yo[i], sigYo[i], 1./sigYo[i]**2,
                                        Phio[i], sigPhio[i], 1./sigPhio[i]**2)
            else:
                self._H[h].add_observation(entering, Frameo[i],
                                        Xo[i], sigXo[i], 1./sigXo[i]**2,
                                        Yo[i], sigYo[i], 1./sigYo[i]**2,
                                        Phio[i], sigPhio[i], 1./sigPhio[i]**2)

    def _spindle_beam_plane_normal(self):
        '''return a unit vector that when placed at the origin of reciprocal
        space, points to the hemisphere of the Ewald sphere
        in which reflections cross from inside to outside of the sphere'''

        # NB vector in +ve Y direction when using imgCIF coordinate frame
        return matrix.col(self._beam.get_s0()).cross(
                        matrix.col(self._gonio.get_rotation_axis())).normalize()

    def _remove_excluded_obs(self, obs_data):
        '''Filter observations that fail certain conditions.

        This includes outlier rejection plus rejection of reflections
        too close to the spindle'''

        # TODO Add outlier rejection. Should use robust statistics.
        # See notes on M-estimators from Garib (when I have them)

        inc = [(h, s, x, sx, y, sy, p, sp) for
               (h, s, x, sx, y, sy, p, sp) in obs_data if
               self._inclusion_test(h, s, self._vecn)]

        return tuple(inc)

    def _inclusion_test(self, h, s, vecn):
        '''Test observation h for inclusion'''

        # At the moment we use a simple catch-all for 'dodgy' observations,
        # removing any that form too small an angle with the spindle_beam plane
        # of the Ewald sphere. This is currently hard-coded. A better test
        # would be adaptive, but to what?

        # reject reflections closer than 6 degrees to the spindle-beam plane
        test = s.accute_angle(vecn) < 1.466077

        return test

    def _create_working_set(self, nref_per_degree):
        '''Make a subset of data for use in refinement'''

        working_data = self._obs_data
        if nref_per_degree:
            temp = self._scan.get_oscillation_range(deg=True)
            width = abs(temp[1] - temp[0])
            sample_size = int(nref_per_degree * width)
            if sample_size < len(working_data):
                self._sample_size = sample_size
                working_data = random.sample(working_data,
                                             self._sample_size)
        return(working_data)

    def get_sample_size(self):
        '''Return the number of observations in the working set to be
        used for refinement'''

        return self._sample_size

    def get_matches(self):
        '''For every observation matched with a prediction return all data'''

        l = [obs for v in self._H.values() for obs in v.obs if obs.is_matched]

        if self._verbosity > 2 and len(l) > 20:
            print "Listing predictions matched with observations for the first 20 reflections:"
            print "H, K, L, Xresid, Yresid, Phiresid"
            fmt = "(%3d, %3d, %3d) %5.3f %5.3f %6.4f"
            for i in xrange(20):
                e = l[i]
                msg = fmt % tuple(e.H + (e.Xresid,
                                 e.Yresid,
                                 e.Phiresid))
                print msg
            print

        return l

    def strip_unmatched_observations(self):
        '''
        Delete observations from the manager that are not matched to a
        prediction. Typically used once, after the first update of
        predictions.
        '''

        for k, v in self._H.items():

            v.remove_unmatched_obs()

            # if no observations left, delete the hkl from the dict
            if len(v) == 0:
                del self._H[k]

        return

    def get_indices(self):
        '''Get the unique indices of all observations in the manager'''

        return flex.miller_index(self._H.keys())

    def get_obs(self, h):
        '''Get the observations of a particular hkl'''

        return self._H[h]

    def get_accepted_reflection_count(self):
        '''Get the number of reflections currently to be used for refinement'''

        return sum(v.get_num_pairs() for v in self._H.values())

    def reset_accepted_reflections(self):
        '''Reset all observations to use=False in preparation for a new set of
        predictions'''

        for v in self._H.values(): v.reset_predictions()
