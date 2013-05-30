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

        # get predictions from the current model and the list of observed hkls

        # FIXME Reorganise this so that there is a loop over all observations
        # that are flagged to be used in refinement. For each observation,
        # predict its rotation angles (what to do for an observation that fails
        # to produce predicted angles?) and the impact positions. Update these
        # each loop cycle. This removes the need for reflection matching,
        # because the observations are explicitly stepped through. Use the
        # DIALS prediction code for this. If it is too slow doing this loop in
        # Python, then move it to C++ (but that will require the observation
        # data structure to be in C++ as well). This needs reworking of the
        # ReflectionManager and other classes too.

        # update the reflection_predictor with current geometry
        self._reflection_predictor.update()

        # predict for all observations in the manager
        predictions = self._reflection_predictor.predict(
                                                self._H.get_indices())

        # obtain the impact positions, currently assuming all
        # reflections intersect panel 0
        impacts = ray_intersection(self._detector, predictions, panel=0)

        # update the ReflectionManager
        self._H.update_predictions(impacts)

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

        self._gradients = self._prediction_parameterisation.get_multi_gradients(self._matches)

        # return residuals and weights as 1d flex.double vectors
        # that is, unroll X, Y and Phi residuals for each match.
        nelem = len(self._matches) * 3
        residuals = flex.double(nelem)
        jacobian_t = flex.double(flex.grid(
            len(self._prediction_parameterisation), nelem))
        weights = flex.double(nelem)

        for i, (m, g) in enumerate(zip(self._matches, self._gradients)):
            residuals[3*i] = m.Xresid
            residuals[3*i + 1] = m.Yresid
            residuals[3*i + 2] = m.Phiresid

            # are these the right weights? Or inverse, or sqrt?
            weights[3*i] = m.weightXo
            weights[3*i + 1] = m.weightYo
            weights[3*i + 2] = m.weightPhio

            #print "X residual, weight = ", residuals[3*i], weights[3*i]
            #print "Y residual, weight = ", residuals[3*i + 1], weights[3*i + 1]
            #print "Phi residual, weight = ", residuals[3*i + 2], weights[3*i + 2]

            dX_dp, dY_dp, dPhi_dp = zip(*g)
            # fill jacobian elements here.
            # g is a nparam length list, each element of which is a triplet of
            # values, (dX/dp_n, dY/dp_n, dPhi/dp_n)

            #print "dX/dp = ", dX_dp
            #print "dY/dp = ", dY_dp
            #print "dPhi/dp = ", dPhi_dp

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

        # for each reflection, get the gradients wrt each parameter
        #self._gradients = [self._prediction_parameterisation.get_gradients(m.H,
        #                        m.Sc, m.Phic) for m in self._matches]
        self._gradients = self._prediction_parameterisation.get_multi_gradients(self._matches)

        for m, grads in zip(self._matches, self._gradients):

            for j, (grad_X, grad_Y, grad_Phi) in enumerate(grads):
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
        for m, grads in zip(self._matches, self._gradients):

            for j, (grad_X, grad_Y, grad_Phi) in enumerate(grads):
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

class ObservationPrediction(object):
    '''A helper class for the reflection manager to contain information about
    the unique observations of particular hkl paired with the currently
    predicted values.

    Reflections are either exiting or entering the Ewald sphere. This is
    calculated by the reflection manager and passed in, and used here to
    identify which observations to update with a new prediction'''

    def __init__(self, H,
                       Xo, sigXo, weightXo,
                       Yo, sigYo, weightYo,
                       Phio, sigPhio, weightPhio, entering):
        '''initialise from one observation'''
        assert isinstance(entering, bool)
        self.H = H
        self.Xo = [Xo]
        self.sigXo = [sigXo]
        self.weightXo = [weightXo]
        self.Yo = [Yo]
        self.sigYo = [sigYo]
        self.weightYo = [weightYo]
        self.Phio = [Phio]
        self.sigPhio = [sigPhio]
        self.weightPhio = [weightPhio]
        self.Xc = [None]
        self.Yc = [None]
        self.Phic = [None]
        self.Sc = [None]
        self.entering = [entering]
        self.use = [False]
        self.num_obs = 1

    def get_num_pairs(self):
        '''Count the number of observations that are paired with a prediction'''
        return sum(1 for x in self.use if x)

    #def reset_predictions(self):
    #    '''Set the use flag to false for all observations, so that after
    #    updating, any observations that still do not have a prediction are
    #    flagged to be removed from calculation of residual and gradients.'''
    #    self.use = [False] * self.num_obs

    def add_observation(self, Xo, sigXo, weightXo,
                              Yo, sigYo, weightYo,
                              Phio, sigPhio, weightPhio, entering):
        '''Add another observation of this reflection'''
        self.Xo.append(Xo)
        self.sigXo.append(sigXo)
        self.weightXo.append(weightXo)
        self.Yo.append(Yo)
        self.sigYo.append(sigYo)
        self.weightYo.append(weightYo)
        self.Phio.append(Phio)
        self.sigPhio.append(sigPhio)
        self.weightPhio.append(weightPhio)
        self.entering.append(entering)
        self.Xc.append(None)
        self.Yc.append(None)
        self.Phic.append(None)
        self.Sc.append(None)
        self.use.append(False)
        self.num_obs += 1

    def update_prediction(self, ref, first_update = False):
        '''Update the observations with a new prediction.'''

        # Current behaviour is to update all observations in the same
        # hemisphere as this prediction. This implies that all reflections
        # in the ReflectionManager come from a model with the same parameter
        # values. This is not the case in reality, as repeat observations of one
        # reflection may be at multiples of 2.*pi from each other, with
        # increasing dose delivered, therefore different cells, etc.
        # Different reflection managers will be needed for each local
        # refinement, where parameter values may differ. In future, it may be
        # better to have one global reflection manager that is able to identify
        # which of multiple observations in one hemisphere a new prediction
        # corresponds to.
        # Anyway, this will have to be addressed to accommodate time dependent
        # parameterisation of models.
        to_update = [ref.entering == n for n in self.entering]

        for i, t in enumerate(to_update):

            if t: # update the prediction for this observation
                Xc, Yc = ref.image_coord_mm
                self.Xc[i] = Xc
                self.Yc[i] = Yc
                # do not wrap around multiples of 2*pi; keep the full rotation
                # from zero to differentiate repeat observations.
                resid = ref.rotation_angle - (self.Phio[i] % TWO_PI)
                self.Phic[i] = self.Phio[i] + resid
                self.Sc[i] = matrix.col(ref.beam_vector)

                #if not self.use[i]:
                # set use flags only for observations that have a prediction
                # on the first cycle
                if first_update:
                    self.use[i] = True

class ObsPredMatch:
    '''A bucket class containing data for a prediction that has been
    matched to an observation'''

    def __init__(self, hkl, Xo, weightXo,
                            Yo, weightYo,
                            Phio, weightPhio,
                            Xc, Yc, Phic,
                            Sc):
        self.H = hkl
        self.Xresid = Xc - Xo
        self.Xresid2 = self.Xresid**2
        #self.Xo = Xo
        self.weightXo = weightXo
        self.Yresid = Yc - Yo
        self.Yresid2 = self.Yresid**2
        #self.Yo = Yo
        self.weightYo = weightYo
        self.Phiresid = Phic - Phio
        self.Phiresid2 = self.Phiresid**2
        #self.Phio = Phio
        self.weightPhio = weightPhio
        #self.Xc = Xc
        #self.Yc = Yc
        self.Phic = Phic
        self.Sc = Sc

class ReflectionManager(object):
    '''A class to maintain information about observed and predicted
    reflections for refinement.'''

    def __init__(self, Ho, So,
                       Xo, sigXo,
                       Yo, sigYo,
                       Phio, sigPhio,
                       beam, gonio, verbosity=0):

        # check the observed values
        Ho = list(Ho)
        So = list(So)
        Xo = list(Xo)
        sigXo = list(sigXo)
        Yo = list(Yo)
        sigYo = list(sigYo)
        Phio = list(Phio)
        sigPhio = list(sigPhio)
        assert(len(So) == \
               len(Xo) == \
               len(sigXo) == \
               len(Yo) == \
               len(sigYo) == \
               len(Phio) == \
               len(sigPhio) == \
               len(Ho))

        # set verbosity
        self._verbosity = verbosity

        # track whether this is the first update of predictions or not
        self._first_update = True

        # keep references to the beam and goniometer models (for reflection
        # exclusion test)
        self._beam = beam
        self._gonio = gonio

        # find vector normal to the spindle-beam plane. In principle this could
        # change during refinement, if we do not fix the beam in a plane.
        # However, for now we determine it once to categorise all observations
        # and predictions into those entering and those leaving the Ewald sphere
        self._vecn = self._spindle_beam_plane_normal()

        # exclude reflections that fail inclusion criteria
        obs_data = zip(Ho, So, Xo, sigXo, Yo, sigYo, Phio, sigPhio)
        (Ho, So, Xo, sigXo, Yo,
         sigYo, Phio, sigPhio) = self._remove_excluded_obs(obs_data)

        # store observation information in a dict of observation-prediction
        # pairs (prediction information will go in here later)
        self._H = {}
        for i, h in enumerate(Ho):
            entering = So[i].dot(self._vecn) < 0.
            if h not in self._H:
                self._H[h] = ObservationPrediction(h,
                                        Xo[i], sigXo[i], 1./sigXo[i]**2,
                                        Yo[i], sigYo[i], 1./sigYo[i]**2,
                                        Phio[i], sigPhio[i], 1./sigPhio[i]**2,
                                        entering)
            else:
                self._H[h].add_observation(Xo[i], sigXo[i], 1./sigXo[i]**2,
                                        Yo[i], sigYo[i], 1./sigYo[i]**2,
                                        Phio[i], sigPhio[i], 1./sigPhio[i]**2,
                                        entering)

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

        return tuple(zip(*inc))

    def _inclusion_test(self, h, s, vecn):
        '''Test observation h for inclusion'''

        # At the moment we use a simple catch-all for 'dodgy' observations,
        # removing any that form too small an angle with the spindle_beam plane
        # of the Ewald sphere. This is currently hard-coded. A better test
        # would be adaptive, but to what?

        # reject reflections closer than 6 degrees to the spindle-beam plane
        test = s.accute_angle(vecn) < 1.466077

        return test

    def get_matches(self):
        '''For every observation matched with a prediction return all data'''

        l = []
        for hkl, v in self._H.items():

            for i, u in enumerate(v.use):

                if u: l.append(ObsPredMatch(hkl, v.Xo[i], v.weightXo[i],
                                            v.Yo[i], v.weightYo[i],
                                            v.Phio[i], v.weightPhio[i],
                                            v.Xc[i], v.Yc[i], v.Phic[i],
                                            v.Sc[i]))

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

    def get_indices(self):
        '''Get the unique indices of all observations in the manager'''

        return flex.miller_index(self._H.keys())

    def get_accepted_reflection_count(self):
        '''Get the number of reflections currently to be used for refinement'''

        return sum(v.get_num_pairs() for v in self._H.values())

    def update_predictions(self, predictions):
        '''Update with the latest values for the predictions.

        Any observations that do not have a prediction are flagged to be
        removed from calculation of residual and gradients.'''

        # Remove all existing predictions (i.e. set use flags to False)
        #for v in self._H.values():
        #    v.reset_predictions()

        # Loop over new predictions, updating matches
        for ref in predictions:

            if ref.miller_index in self._H: # found an observation for this prediction

                h = ref.miller_index
                s = matrix.col(ref.beam_vector)
                x, y = ref.image_coord_mm
                p = ref.rotation_angle

                # exclude reflections that fail inclusion criteria
                if not self._inclusion_test(
                        h, s, self._vecn): continue

                self._H[h].update_prediction(ref,
                    first_update = self._first_update)
        self._first_update = False

        return
