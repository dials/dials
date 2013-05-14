#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

from __future__ import division
from dials.algorithms.refinement.parameterisation.model_parameters import *
from math import exp
from scitbx import matrix
from dials.algorithms.refinement import dR_from_axis_and_angle

class ScanVaryingParameterSet(Parameter):
    '''Testing a class for a scan-varying parameter, in which values at rotation
    angle phi may be derived using smoothed interpolation between checkpoint
    values stored here. Externally, this is presented as a set of parameters.

    num_samples is the number of checkpoints. Other arguments are as Parameter.
    '''

    def __init__(self, value, num_samples = 5, axis = None, ptype = None, name = "ScanVaryingParameterSet"):

        Parameter.__init__(self, value, axis, ptype, name)

        assert num_samples >= 2 # otherwise use scan-independent parameterisation
        self._num_samples = num_samples
        self._value = [value] * num_samples
        self._esd = [None] * num_samples
        self._axis = axis
        self._ptype = ptype
        name_stem = [name] * num_samples
        self._name = [e + "_sample%d" % i for i, e in enumerate(name_stem)]
        self._fixed = False

        return

    def __len__(self):
        return self._num_samples

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, val):
        assert len(val) == len(self)
        self._value = val
        self._esd = [None] * len(self)

    @property
    def axis(self):
        return self._axis


class GaussianSmoother(object):
    '''A Gaussian smoother for ScanVaryingModelParameterisations'''

    # Based largely on class SmoothedValue from Aimless.

    # Construct from range of raw unnormalised coordinate & number of sample intervals
    # Set smoothing values to defaults, Nav = 3
    def __init__(self, x_range, num_intervals):

        self._x0 = x_range[0] # coordinate of z = 0
        self._nsample = num_intervals # number of intervals
        assert self._nsample > 0
        if self._nsample == 1:
            self._nvalues = 2
        elif self._nsample == 2:
            self._nvalues = 3
        else:
            self._nvalues = self._nsample + 2

        # smoothing spacing
        self._spacing = (x_range[1] - x_range[0]) / float(self._nsample)

        # the values are actually held by ScanVaryingParameterSet classes, but
        # we need the positions
        if self._nvalues == 2:
            self._positions = [1.0, 2.0]
        if self._nvalues == 3:
            self._positions = [0.0, 1.0, 2.0]
        else:
            self._positions = [e - 0.5 for e in range(self._nvalues)]

        # set default smoothing parameters
        self.set_smoothing(3, -1.0)

    def set_smoothing(self, num_average, sigma):

        '''Set smoothing values:

        naverage: number of points included in each calculation
        sigma: width of the Gaussian used for smoothing.

        If sigma < 0, set to "optimum" (!) (or at least suitable) value from
        num_average '''

        self._naverage = num_average
        if self._naverage > self._nvalues:
            self._naverage = self._nvalues
        self._half_naverage = self._naverage / 2.0
        self._sigma = sigma

        if self._naverage < 1 or self._naverage > 5:
            msg = "num_average must be between 1 & 5"
            raise ValueError, msg

        if sigma < 0.0:
            #Default values 0.65, 0.7, 0.75, 0.8 for nav = 2,3,4,5
            self._sigma = 0.65 + 0.05 * (self._naverage - 2)

    # Return number of values
    def num_values(self):
        return self._nvalues

    # Return number of sample intervals
    def num_samples(self):
        return self._nsample

    # Return interpolated value of param at point, original unnormalised
    # coordinate. Also return the weights at each position.
    def value_weight(self, x, param):
        pass

        weight = [0.0] * len(self._positions)

        # normalised coordinate
        z = (x - self._x0) / self._spacing
        sumwv = 0.0
        sumweight = 0.0

        # get values
        values = param.value

        if self._nvalues <= 3:
            i1 = 0
            i2 = self._nvalues
        else: # in this case, 1st point in array (index 0) is at position -0.5
            # find the nearest naverage points that bracket z
            i1 = int(round(z - self._half_naverage)) + 1
            i2 = i1 + self._naverage
            if i1 < 0: # beginning of range
                i1 = 0
                i2 = max(2, i2) # ensure a separation of at least 2
            if i2 > self._nvalues:
                i2 = self._nvalues
                i1 = min(i1, self._nvalues - 2) # ensure a separation of at least 2

        # now do stuff
        for i in range(i1, i2):

            ds = (z - self._positions[i]) / self._sigma
            weight[i] = exp(-ds*ds)
            sumwv += weight[i] * values[i]
            sumweight  += weight[i]

        if sumweight > 0.0:
            value = sumwv / sumweight;
        else:
            value = 0

        return value, weight, sumweight

    # Return number of points averaged
    def num_average(self):
        return self._naverage

    # Return sigma smoothing factor
    def sigma(self):
        return self._sigma

    # Return spacing
    def spacing(self):
        return self._spacing

    # Return positions
    def positions(self):
        return self._positions



class ScanVaryingModelParameterisation(ModelParameterisation):
    '''Extending ModelParameterisation to deal with ScanVaryingParameterSets.

    For simplicity at this stage it is decreed that a
    ScanVaryingModelParameterisation consists only of ScanVaryingParameterSets.
    There is no combination with normal Parameters. This could be changed later,
    but there may be no reason to do so, hence starting with this simpler
    design'''

    # The initial state is here equivalent to the initial state of the
    # time static version of the parameterisation, as it is assumed that we
    # start with a flat model wrt rotation angle.

    def __init__(self, models, initial_state, param_sets, smoother):
        assert(isinstance(param_sets, list))
        self._initial_state = initial_state
        self._models = models
        self._param_sets = list(param_sets)
        self._num_sets = len(self._param_sets)
        self._set_len = len(param_sets[0])
        self._total_len = self._set_len * self._num_sets

        # ensure all internal parameter sets have the same number of parameters
        for param in self._param_sets[1:]: assert len(param) == self._set_len

        #self._dstate_dp = [None] * self._total_len
        #self._pset_fixed = [False] * self._num_sets

        # Link up with an object that will perform the smoothing.
        self._smoother = smoother
        assert self._smoother.num_values() == self._set_len

        return

    #def __len__(self):
    #    return len(self._param_sets)

    def num_free(self):
        '''the number of free parameters'''
        return sum(not x.get_fixed() for x in self._param_sets) * self._set_len

    def num_total(self):
        '''the total number of parameters, both fixed and free'''
        return self._total_len

    def compose(self):
        '''compose the current model state from its initial state and its
        parameter list. Also calculate the derivatives of the state wrt
        each parameter in the list. Should be called automatically once
        parameters are updated, e.g. at the end of each refinement cycle'''

        raise RuntimeError, 'implement me'

    def get_p(self, only_free = True):
        '''export the values of the internal list of parameters as a
        sequence of floats.

        If only_free, the values of fixed parameters are filtered from the
        returned list. Otherwise all parameter values are returned'''

        if only_free:
            return [x for e in self._param_sets \
                    if not e.get_fixed() for x in e.value]
            #return [x for e, f in zip(self._param_sets, self._pset_fixed) \
            #        if not f for x in e.value]

        else:
            return [x for e in self._param_sets for x in e.value]

    def get_pnames(self, only_free = True):
        '''export the names of the internal list of parameters

        If only_free, the names of fixed parameters are filtered from the
        returned list. Otherwise all parameter names are returned'''

        # FIXME combine functionality with get_p by returning a named, ordered
        # list
        if only_free:
            return [x for e in self._param_sets \
                    if not e.get_fixed() for x in e.name]
            #return [x for e, f in zip(self._param_sets, self._pset_fixed) \
            #        if not f for x in e.name]

        else:
            return [x for e in self._param_sets for x in e.name]

    def set_p(self, vals):
        '''set the values of the internal list of parameters from a
        sequence of floats.

        First break the sequence into sub sequences of the same length
        as the _set_len.

        Only free parameter sets can have values assigned, therefore the
        length of vals must equal the value of num_free'''

        assert(len(vals) == self.num_free())
        i = 0
        for p in self._param_sets:
            if not p.get_fixed(): # only set the free parameter sets
                new_vals = vals[i:i+self._set_len]
                p.value = new_vals
                i += self._set_len

        # compose with the new parameter values
        #self.compose()

        return

    def get_fixed(self):
        '''return the list determining whether each parameter set is fixed or not'''
        return [p.get_fixed() for p in self._param_sets]

    def set_fixed(self, fix):
        '''set parameter sets to be fixed or free from a list'''

        assert(len(fix) == len(self._param_sets))

        for f, p in zip(fix, self._param_sets):
            if f: p.fix()
            else: p.unfix()

    def get_state(self, t):
        '''return the current state of the model under parameterisation
        at time t. This is required, for example, by the calculation
        of finite difference gradients.'''

        # To be implemented by the derived class, where it is clear what aspect
        # of the model under parameterisation is considered its state. The
        # type of this result should match the type of one element of the return
        # value of get_ds_dp.
        raise RuntimeError, 'implement me'

    def get_ds_dp(self, t, only_free = True):
        '''get a list of derivatives of the state wrt each parameter, as
        a list in the same order as the internal list of parameters. Evaluate
        each scan-dependent parameter at coordinate t, corresponding to
        the original, unnormalised coordinates used to set up the smoother
        (t will most likely be along the dimension of image number).

        If only_free, the derivatives with respect to fixed parameters are
        omitted from the returned list. Otherwise a list for all parameters is
        returned, with values of 0.0 for the fixed parameters'''

        # To be implemented by the derived class. A list _dstate_dp cannot
        # reasonably be stored here, because the derivatives are time
        # dependent, so must be calculated for each requested point t
        raise RuntimeError, 'implement me'

        # make a list that tells whether each parameter is fixed
#        p_fixed = [x for f in self.get_fixed() for x in [f] * self._set_len]
#
#        if only_free:
#
#            return [e for e, f in zip(self._dstate_dp,
#                                      p_fixed) if not f]
#
#        else:
#            return [0. * e if f else e for e, f in zip(self._dstate_dp,
#                                                       p_fixed)]


class ScanVaryingCrystalOrientationParameterisation(ScanVaryingModelParameterisation):
    '''A work-in-progress time-dependent parameterisation for crystal
    orientation, with angles expressed in mrad'''

    def __init__(self, crystal, t_range, num_intervals):

        # The state of a scan varying crystal orientation parameterisation
        # is an orientation
        # matrix '[U](t)', expressed as a function of 'time' t (which could
        # actually be measured by image number in a sequential scan)

        # The initial state is a snapshot of the crystal orientation
        # at the time of initialisation '[U0]', which independent of time.
        # Future states are composed by
        # rotations around axes of the phi-axis frame by Tait-Bryan angles.
        #
        # [U](t) = [Phi3](t)[Phi2](t)[Phi1](t)[U0]

        # Set up the smoother
        smoother = GaussianSmoother(t_range, num_intervals)
        nv = smoother.num_values()

        ### Set up the initial state
        istate = crystal.get_U()

        ### Set up the parameters
        phi1 = ScanVaryingParameterSet(0.0, nv,
                            matrix.col((1., 0., 0.)), 'angle', 'Phi1')
        phi2 = ScanVaryingParameterSet(0.0, nv,
                            matrix.col((0., 1., 0.)), 'angle', 'Phi2')
        phi3 = ScanVaryingParameterSet(0.0, nv,
                            matrix.col((1., 0., 0.)), 'angle', 'Phi3')

        # build the list of parameter sets in a specific, maintained order
        p_list = [phi1, phi2, phi3]

        # set up the list of model objects being parameterised (here
        # just a single crystal model)
        models = [crystal]

        # set up the base class
        ScanVaryingModelParameterisation.__init__(self, models, istate,
                                                  p_list, smoother)

        # call compose to calculate all the derivatives
        #FIXME cannot use compose to calculate derivatives as they are
        # time dependent now
        #self.compose()

    def get_ds_dp(self, t, only_free = True):
        '''calculate derivatives for model at time t'''

        # Extract orientation from the initial state
        U0 = self._initial_state

        # extract parameter sets from the internal list
        phi1_set, phi2_set, phi3_set = self._param_sets

        # extract angles and other data at time t using the smoother
        phi1, phi1_weights, phi1_sumweights = self._smoother.value_weight(t, phi1_set)
        phi2, phi2_weights, phi2_sumweights = self._smoother.value_weight(t, phi2_set)
        phi3, phi3_weights, phi3_sumweights = self._smoother.value_weight(t, phi3_set)

        # calculate derivatives of angles wrt underlying parameters.
        # FIXME write up notes in orange notebook
        dphi1_dp = [e / phi1_sumweights for e in phi1_weights]
        dphi2_dp = [e / phi2_sumweights for e in phi2_weights]
        dphi3_dp = [e / phi3_sumweights for e in phi3_weights]

        # convert angles to radians
        phi1rad, phi2rad, phi3rad = (phi1 / 1000., phi2 / 1000.,
                                     phi3 / 1000.)

        # compose rotation matrices and their first order derivatives wrt angle
        Phi1 = (phi1_set.axis).axis_and_angle_as_r3_rotation_matrix(phi1rad, deg=False)
        dPhi1_dphi1 = dR_from_axis_and_angle(phi1_set.axis, phi1rad, deg=False)

        Phi2 = (phi2_set.axis).axis_and_angle_as_r3_rotation_matrix(phi2rad, deg=False)
        dPhi2_dphi2 = dR_from_axis_and_angle(phi2_set.axis, phi2rad, deg=False)

        Phi3 = (phi3_set.axis).axis_and_angle_as_r3_rotation_matrix(phi3rad, deg=False)
        dPhi3_dphi3 = dR_from_axis_and_angle(phi3_set.axis, phi3rad, deg=False)

        Phi21 = Phi2 * Phi1
        Phi321 = Phi3 * Phi21

        ### Compose new state

        #newU = Phi321 * U0
        #self._models[0].set_U(newU)

        ### calculate derivatives of the state wrt angle, convert back to mrad
        dU_dphi1 = Phi3 * Phi2 * dPhi1_dphi1 * U0 / 1000.
        dU_dphi2 = Phi3 * dPhi2_dphi2 * Phi1 * U0 / 1000.
        dU_dphi3 = dPhi3_dphi3 * Phi21 * U0 / 1000.

        # calculate derivatives of state wrt underlying parameters
        dU_dp1 = [dU_dphi1 * e for e in dphi1_dp]
        dU_dp2 = [dU_dphi2 * e for e in dphi2_dp]
        dU_dp3 = [dU_dphi3 * e for e in dphi3_dp]

        # return concatenated list of derivatives
        return dU_dp1 + dU_dp2 + dU_dp3

    def get_state(self, t):

        '''Return crystal orientation matrix [U] at time t'''

        # Extract orientation from the initial state
        U0 = self._initial_state

        # extract parameter sets from the internal list
        phi1_set, phi2_set, phi3_set = self._param_sets

        # extract angles and other data at time t using the smoother
        phi1, phi1_weights, phi1_sumweights = self._smoother.value_weight(t, phi1_set)
        phi2, phi2_weights, phi2_sumweights = self._smoother.value_weight(t, phi2_set)
        phi3, phi3_weights, phi3_sumweights = self._smoother.value_weight(t, phi3_set)

        # convert angles to radians
        phi1rad, phi2rad, phi3rad = (phi1 / 1000., phi2 / 1000.,
                                     phi3 / 1000.)

        # compose rotation matrices
        Phi1 = (phi1_set.axis).axis_and_angle_as_r3_rotation_matrix(phi1rad, deg=False)
        Phi2 = (phi2_set.axis).axis_and_angle_as_r3_rotation_matrix(phi2rad, deg=False)
        Phi3 = (phi3_set.axis).axis_and_angle_as_r3_rotation_matrix(phi3rad, deg=False)

        Phi21 = Phi2 * Phi1
        Phi321 = Phi3 * Phi21

        ### Compose new state

        newU = Phi321 * U0
        #self._models[0].set_U(newU)
        # get U(t)
        return newU
