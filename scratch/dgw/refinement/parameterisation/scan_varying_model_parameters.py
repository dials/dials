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

class ScanVaryingParameter(Parameter):
    '''Testing a class for a scan-varying parameter, in which values at rotation
    angle phi may be derived using smoothed interpolation between checkpoint
    values stored here.

    num_samples is the number of checkpoints. Other arguments are as Parameter.
    '''

    def __init__(self, value, num_samples = 5, axis = None, ptype = None, name = "ScanVaryingParameter"):

        Parameter.__init__(self, value, axis, ptype, name)

        assert num_samples >= 5 # could be lower
        self._num_samples = num_samples
        self._value = [value] * num_samples
        self._esd = [None] * num_samples
        self._axis = axis
        self._ptype = ptype
        name_stem = [name] * num_samples
        self._name = [e + "_sample%d" % i for i, e in enumerate(name_stem)]

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
    def __init__(self, phi_range, num_intervals):

        self._x0 = phi_range[0] # coordinate of z = 0
        self._nsample = num_intervals # number of intervals
        assert self._nsample > 0 # otherwise use scan-independent parameterisation
        if self._nsample == 1:
            self._nvalues = 2
        elif self._nsample == 2:
            self._nvalues = 3
        else:
            self._nvalues = self._nsample + 2

        # smoothing spacing
        self._spacing = (phi_range[1] - phi_range[0]) / float(self._nsample)

        # the values are actually held by ScanVaryingParameter classes, but
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

        return value, weight

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
    '''Extending ModelParameterisation to deal with ScanVaryingParameters.

    For simplicity at this stage it is decreed that a
    ScanVaryingModelParameterisation consists only of ScanVaryingParameters.
    There is no combination with normal Parameters. This could be changed later,
    but there may be no reason to do so, hence starting with this simpler
    design'''

    # The initial state is here equivalent to the initial state of the
    # time static version of the parameterisation, as it is assumed that we
    # start with a flat model wrt rotation angle.

    def __init__(self, models, initial_state, param_list, num_samples, smoother):
        assert(isinstance(param_list, list))
        self._initial_state = initial_state
        self._models = models
        self._plist = list(param_list)
        self._dstate_dp = [None] * len(param_list)
        self._pfixed = [False] * len(param_list)

        # Choose the number of checkpoints for each parameter. This could be
        # determined automatically in the scope that creates this object,
        # or via user preferences. It will affect the smoothing properties of
        # the parameterisation.
        self._num_samples = num_samples

        # Link up with an object that will perform the smoothing
        self._smoother = smoother

        return

    #def __len__(self):
    #    return len(self._plist)

    def num_free(self):
        '''the number of free parameters'''
        return len([x for x in self._pfixed if not x])

    def num_total(self):
        '''the total number of parameters, both fixed and free'''
        return len(self._plist)

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
            return [x.value for x, f in zip(self._plist, self._pfixed) if not f]

        else:
            return [x.value for x in self._plist]

    def get_pnames(self, only_free = True):
        '''export the names of the internal list of parameters

        If only_free, the names of fixed parameters are filtered from the
        returned list. Otherwise all parameter names are returned'''

        # FIXME combine functionality with get_p by returning a named, ordered
        # list
        if only_free:
            return [x.name for x, f in zip(self._plist, self._pfixed) if not f]

        else:
            return [x.name for x in self._plist]

    def set_p(self, vals):
        '''set the values of the internal list of parameters from a
        sequence of floats.

        Only free parameters can be set, therefore the length of vals must equal
        the value of num_free'''

        assert(len(vals) == self.num_free())
        for par, new_val in zip((p for p, f in zip(self._plist, self._pfixed) if not f), vals):
            par.value = new_val

        # compose with the new parameter values
        self.compose()

        return

    def get_fixed(self):
        '''return the list determining whether each parameter is fixed or not'''
        return list(self._pfixed)

    def set_fixed(self, fix):
        '''set the list determining whether each parameter is fixed or not'''

        assert(len(fix) == len(self._plist))
        self._pfixed = [True if e else False for e in fix]

    def get_state(self):
        '''return the current state of the model under parameterisation. This is
        required, for example, by the calculation of finite difference
        gradients.'''

        # To be implemented by the derived class, where it is clear what aspect
        # of the model under parameterisation is considered its state. The
        # type of this result should match the type of one element of the return
        # value of get_ds_dp.
        raise RuntimeError, 'implement me'

    def get_ds_dp(self, at_phi, only_free = True):
        '''get a list of derivatives of the state wrt each parameter, as
        a list in the same order as the internal list of parameters. Evaluate
        each scan-dependent parameter at angle at_phi

        If only_free, the derivatives with respect to fixed parameters are
        omitted from the returned list. Otherwise a list for all parameters is
        returned, with values of 0.0 for the fixed parameters'''

        if only_free:
            return [e for e, f in zip(self._dstate_dp,
                                      self._pfixed) if not f]

        else:
            return [0. * e if f else e for e, f in zip(self._dstate_dp,
                                                       self._pfixed)]
