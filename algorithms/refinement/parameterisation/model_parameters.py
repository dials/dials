#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

from __future__ import division

class Parameter(object):
    '''A class to help formalise what a parameter is. A Parameter must
    have a numerical value (either a length or an angle). It may also
    have a vector axis which provides context for what that number
    means.

    Together the values and axes of a set of parameters' can be
    used to compose the state of a model. For example, the value might be
    a rotation angle, with the axis of rotation providing the context.

    A slot is also provided for the estimated standard deviation of the
    value, which may be of use in future. Currently, whenever the
    parameter value is set, the esd is reset to None. So this must be
    set separately, and after the parameter value if it is required'''

    def __init__(self, value, axis = None, ptype = None, name = "Parameter"):
        self._value = value
        self._esd = None
        self._axis = axis
        #assert ptype in ['length', 'angle']
        self._ptype = ptype
        self._name = name
        self._fixed = False

        return

    @property
    def value(self):
        return self._value

    @property
    def name(self):
        return self._name

    @value.setter
    def value(self, val):
        self._value = val
        self._esd = None

    @property
    def esd(self):
        return self._esd

    @esd.setter
    def esd(self, esd):
        self._esd = esd

    @property
    def axis(self):
        return self._axis

    def get_fixed(self):
        return self._fixed

    def fix(self):
        self._fixed = True

    def unfix(self):
        self._fixed = False



class ModelParameterisation(object):
    '''An abstract interface that model elements, such as the detector
    model, the source model, etc. should adhere to in order to compose
    their state from their parameters, access their parameters, and
    derivatives of their state wrt their parameters, taking into account whether
    particular parameters are fixed or free.'''

    def __init__(self, models, initial_state, param_list):
        assert(isinstance(param_list, list))
        self._initial_state = initial_state
        self._models = models
        self._param = list(param_list)
        self._total_len = len(self._param)
        self._dstate_dp = [None] * len(param_list)

        return

    #def __len__(self):
    #    return len(self._param)

    def num_free(self):
        '''the number of free parameters'''

        return sum(not x.get_fixed() for x in self._param)

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

            return [x.value for x in self._param if not x.get_fixed()]

        else:
            return [x.value for x in self._param]

    def get_pnames(self, only_free = True):
        '''export the names of the internal list of parameters

        If only_free, the names of fixed parameters are filtered from the
        returned list. Otherwise all parameter names are returned'''

        # FIXME combine functionality with get_p by returning a named, ordered
        # list

        if only_free:

            return [x.name for x in self._param if not x.get_fixed()]

        else:
            return [x.name for x in self._param]

    def set_p(self, vals):
        '''set the values of the internal list of parameters from a
        sequence of floats.

        Only free parameters can be set, therefore the length of vals must equal
        the value of num_free'''

        assert(len(vals) == self.num_free())

        v = iter(vals)
        for p in self._param:
            if not p.get_fixed(): # only set the free parameters
                p.value = v.next()

        # compose with the new parameter values
        self.compose()

        return

    def get_fixed(self):
        '''return the list determining whether each parameter is fixed or not'''

        return [p.get_fixed() for p in self._param]


    def set_fixed(self, fix):
        '''set parameters to be fixed or free'''

        assert(len(fix) == len(self._param))

        for f, p in zip(fix, self._param):
            if f: p.fix()
            else: p.unfix()

        return

    def get_state(self):
        '''return the current state of the model under parameterisation. This is
        required, for example, by the calculation of finite difference
        gradients.'''

        # To be implemented by the derived class, where it is clear what aspect
        # of the model under parameterisation is considered its state. The
        # type of this result should match the type of one element of the return
        # value of get_ds_dp.
        raise RuntimeError, 'implement me'

    def get_ds_dp(self, only_free = True):
        '''get a list of derivatives of the state wrt each parameter, as
        a list in the same order as the internal list of parameters.

        If only_free, the derivatives with respect to fixed parameters are
        omitted from the returned list. Otherwise a list for all parameters is
        returned, with values of 0.0 for the fixed parameters'''

        if only_free:
            return [ds_dp for ds_dp, p in zip(self._dstate_dp, self._param) \
                    if not p.get_fixed()]

        else:
            return [0. * ds_dp if p.get_fixed() else ds_dp \
                        for ds_dp, p in zip(self._dstate_dp, self._param)]
