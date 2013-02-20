# Copyright (C) (2012) David Waterman, STFC Rutherford Appleton Laboratory, UK.
# This code is developed as part of the DIALS project and is provided for
# testing purposes only

from __future__ import division
from copy import deepcopy

class parameter(object):
    '''A class to help formalise what a parameter is. A parameter must
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

    def __init__(self, value, axis = None, ptype = None):
        self._value = value
        self._esd = None
        self._axis = axis
        #assert ptype in ['length', 'angle']
        self._ptype = ptype

        return

    @property
    def value(self):
        return self._value

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



class model_parameterisation(object):
    '''An abstract interface that model elements, such as the detector
    model, the source model, etc. should adhere to in order to compose
    their state from their parameters, access their parameters, and
    derivatives of their state wrt their parameters, taking into account whether
    particular parameters are fixed or free.'''

    def __init__(self, models, initial_state, param_list):
        assert(isinstance(param_list, list))
        self._initial_state = initial_state
        self._models = models
        self._plist = list(param_list)
        self._dstate_dp = [None] * len(param_list)
        self._pfixed = [False] * len(param_list)

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

    def get_ds_dp(self, only_free = True):
        '''get a list of derivatives of the state wrt each parameter, as
        a list in the same order as the internal list of parameters.

        If only_free, the derivatives with respect to fixed parameters are
        omitted from the returned list. Otherwise a list for all parameters is
        returned, with values of 0.0 for the fixed parameters'''

        if only_free:
            return [e for e, f in zip(self._dstate_dp,
                                      self._pfixed) if not f]

        else:
            return [0. * e if f else e for e, f in zip(self._dstate_dp,
                                                       self._pfixed)]
