#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

from __future__ import absolute_import, division

from scitbx.array_family import flex
import abc

class Parameter(object):
  """A class to help formalise what a parameter is.

  A Parameter must have a numerical value (either a length or an angle). It may
  also have a vector axis which provides context for what that number means.

  Together the values and axes of a set of parameters can be used to compose
  the state of a model. For example, the value might be a rotation angle, with
  the axis of rotation providing the context.

  A slot is also provided for storage of the estimated standard deviation of
  the value, which may be of use in future. Whenever the parameter
  value is set, the esd is reset to None. So this must be set separately, and
  after the parameter value if it is required.
  """

  def __init__(self, value, axis = None, ptype = None, name = "Parameter"):
    """Initialisation for the Parameter class.

    Args:
        value (float): The initial parameter value.
        axis: A 3-element vector giving the axis of action of the parameter.
            Defaults to None.
        ptype (str): A string defining the parameter type.
        name (str): A string defining the parameter name.
    """
    self._value = value
    self._esd = None
    self._axis = axis
    self._ptype = ptype
    self._name = name
    self._fixed = False

    return

  @property
  def value(self):
    """The floating point parameter value."""
    return self._value

  @value.setter
  def value(self, val):
    self._value = val
    self._esd = None

  @property
  def name(self):
    """A string identifying the parameter name."""
    return self._name

  @property
  def esd(self):
    """The floating point parameter estimated standard deviation value."""
    return self._esd

  @esd.setter
  def esd(self, esd):
    self._esd = esd

  @property
  def param_type(self):
    """A string identifying the parameter type."""
    return self._ptype

  @property
  def axis(self):
    """A 3D column vector giving the direction of action of the parameter."""
    return self._axis

  def get_fixed(self):
    """Return a bool describing whether the parameter is fixed or not."""
    return self._fixed

  def fix(self):
    """Fix the parameter"""
    self._fixed = True

  def unfix(self):
    """Unfix (free) the parameter"""
    self._fixed = False

  def __str__(self):
    """Provide a descriptive string representation of the Parameter class"""

    msg = "Parameter " + self.name + ":\n"
    try:
      msg += "    Type: " + self.param_type + "\n"
    except TypeError:
      msg += "    Type: " + str(self.param_type) + "\n"
    try:
      msg += "    Axis: (%5.3f, %5.3f, %5.3f)" % tuple(self.axis) + "\n"
    except TypeError:
      msg += "    Axis: " + str(self.axis) + "\n"
    msg += "    Value: %5.3f" % self.value + "\n"
    try:
      msg += "    Sigma: %5.3f" % self.esd + "\n"
    except TypeError:
      msg += "    Sigma: " + str(self.esd) + "\n"

    return msg

class ModelParameterisation(object):
  """An abstract interface for the parameterisation of a model.

  Parameterisation of experimental objects, such as the detector, the beam,
  etc. should adhere to this interface in order to compose their state from
  their parameters, access their parameters, and derivatives of their state wrt
  their parameters, taking into account whether particular parameters are fixed
  or free.

  It is possible to parameterise a model with multiple states. The first such
  example is a detector with multiple panels. Each panel has its own matrix
  describing its geometrical 'state'. One set of parameters is used to compose
  all states and calculate all derivatives of these states.
  """

  __metaclass__  = abc.ABCMeta

  def __init__(self, model, initial_state, param_list, experiment_ids,
               is_multi_state=False):
    """Initialisation for the ModelParameterisation abstract base class.

    Args:
        model: The experimental model to be parameterised.
        initial_state: The initial value of the state of interest of the model.
        param_list (list): A list of Parameter objects.
        experiment_ids (list): A list of integer experiment IDs that this
            parameterisation affects.
        is_multi_state (bool): A flag to indicate whether the parameterisation
            is for a multi-state model.
    """
    assert(isinstance(param_list, list))
    self._initial_state = initial_state
    self._model = model
    self._param = list(param_list)
    self._total_len = len(self._param)
    self._num_free = None
    self._dstate_dp = [None] * len(param_list)
    self._is_multi_state = is_multi_state
    if is_multi_state:
      self._multi_state_derivatives = [None] * len(model)
    self._exp_ids = experiment_ids
    self._null_state = self.get_state() * 0.0
    return

  def is_multi_state(self):
    """Query whether this is a multi-state parameterisation or not (e.g. True
    for a multi-panel detector parameterisation)"""

    return self._is_multi_state

  def num_free(self):
    """Get the number of free parameters"""

    if self._num_free is None:
      self._num_free = sum(not x.get_fixed() for x in self._param)
    return self._num_free

  def num_total(self):
    """Get the total number of parameters, both fixed and free"""
    return self._total_len

  def get_experiment_ids(self):
    """Get the list of experiment IDs of parameterised experiments."""
    return self._exp_ids

  def get_model(self):
    """Get the model that is parameterised"""
    return self._model

  @abc.abstractmethod
  def compose(self):
    """Compose the current model state.

    Using the initial state and the current parameter list, compose the current
    state of the model. Also calculate the derivatives of the state wrt each
    parameter in the list. Intended to be called automatically once parameters
    are updated, e.g. at the end of each refinement cycle.
    """

    # Specific methods for each model are defined by derived classes.
    pass

  def get_params(self, only_free = True):
    """Get the internal list of parameters. It is intended that this function
    be used for reporting parameter attributes, not for modifying them.

    Args:
        only_free (bool): Whether fixed parameters should be filtered from the
            returned list, or all parameters returned.
    """

    if only_free:

      return [x for x in self._param if not x.get_fixed()]

    else:
      return [x for x in self._param]

  def get_param_vals(self, only_free = True):
    """Get the values of the internal list of parameters as a sequence of
    floats.

    Args:
        only_free (bool): Whether fixed parameters should have their values
        filtered from the returned list, or all parameter values returned.
    """

    if only_free:

      return [x.value for x in self._param if not x.get_fixed()]

    else:
      return [x.value for x in self._param]

  def get_param_names(self, only_free = True):
    """Get the list of names of the parameters.

    Args:
        only_free (bool): Whether fixed parameters should have their names
        filtered from the returned list, or all parameter names returned.
    """

    if only_free:
      return [x.name for x in self._param if not x.get_fixed()]

    else:
      return [x.name for x in self._param]

  def set_param_vals(self, vals):
    """Set the values of the internal list of parameters from a sequence of
    floats.

    Args:
        vals (list): A list of floating point parameter values, equal in length
            to the number of free parameters.
    """

    assert(len(vals) == self.num_free())

    v = iter(vals)
    for p in self._param:
      if not p.get_fixed(): # only set the free parameters
        p.value = next(v)
        p.esd = None

    # compose with the new parameter values
    self.compose()

    return

  def set_param_esds(self, esds):
    """Set the estimated standard deviations of the internal list of parameters
    from a sequence of floats.

    Args:
        esds (list): A list of floating point parameter esd values, equal in
            length to the number of free parameters.
    """

    assert(len(esds) == self.num_free())

    v = iter(esds)
    for p in self._param:
      if not p.get_fixed(): # only set the free parameters
        p.esd = next(v)

    return

  def get_fixed(self):
    """Get the list determining whether each parameter is fixed or not"""

    return [p.get_fixed() for p in self._param]


  def set_fixed(self, fix):
    """Set parameters to be fixed or free.

    Args:
        fix (list): A list of bools, determining whether each parameter is
            fixed or free.
    """

    assert(len(fix) == len(self._param))

    for f, p in zip(fix, self._param):
      if f: p.fix()
      else: p.unfix()

    # reset the cached number of free parameters
    self._num_free = None

    return

  @abc.abstractmethod
  def get_state(self, multi_state_elt=None):
    """Get the current state of the model under parameterisation.

    This is required, for example, by the calculation of finite
    difference gradients.

    Args:
        multi_state_elt (bool): For a multi-state parameterisation, the
            requested state is selected by an integer array index. Defaults to
            None.
    """

    # To be implemented by the derived class, where it is clear what aspect
    # of the model under parameterisation is considered its state. The
    # type of this result should match the type of one element of the return
    # value of get_ds_dp.
    pass

  def get_ds_dp(self, only_free = True, multi_state_elt=None, use_none_as_null=False):
    """Get a list of derivatives of the state wrt each parameter.

    This is returned as a list in the same order as the internal list of
    parameters.

    The internal list of derivatives may use None for null elements. By default
    these are converted to the null state, but optionally these may remain None
    to detect them easier and avoid doing calculations on null elements.

    Args:
        only_free (bool): whether the derivatives with respect to fixed
            parameters are omitted from the returned list. Otherwise a list for
            all parameters is returned, with null values for the fixed
            parameters.
        multi_state_elt (bool): For a multi-state parameterisation, the
            requested derivative of the state is selected by an integer array
            index. Defaults to None.
        use_none_as_null (bool): See method description. Defaults to False.
    """

    if use_none_as_null:
      null = None
    else:
      null = self._null_state

    if self._is_multi_state:
      if multi_state_elt is None:
        raise ValueError('Must provide multi_state_elt for a multi-state parameterisation')
      else:
        # copy results for the state of interest
        self._dstate_dp = self._multi_state_derivatives[multi_state_elt]

    if only_free:
      grads = [ds_dp for ds_dp, p in zip(self._dstate_dp, self._param) if not p.get_fixed()]
    else:
      grads = self._dstate_dp

    return [null if e is None else e for e in grads]


  def calculate_state_uncertainties(self, var_cov):
    """Given a variance-covariance array for the parameters of this model,
    propagate those estimated errors into the uncertainties of the model state

    Args:
        var_cov: A variance-covariance matrix for the parameters of this model.
    """

    grads = []
    if self._is_multi_state:
      i = 0
      while True:
        try:
          grads.append(self.get_ds_dp(multi_state_elt=i))
          i += 1
        except IndexError:
          break
    else:
      grads.append(self.get_ds_dp())

    if len(grads[0]) == 0: return

    # the jacobian is the m*n matrix of partial derivatives of the m state
    # elements wrt the n parameters
    state_covs = []
    for grads_one_state in grads:
      jacobian_t = flex.double(e for g in grads_one_state for e in g)
      jacobian_t.reshape(flex.grid(len(grads_one_state),
                                   len(grads_one_state[0].elems)))

      # propagation of errors takes the variance-covariance matrix of parameters,
      # along with the jacobian mapping changes in parameter values to changes
      # in the model state elements, to calculate an approximate variance-
      # covariance matrix of the state elements. That is, state_cov is the
      # matrix product: jacobian * var_cov * jacobian_t
      tmp = var_cov.matrix_multiply(jacobian_t)
      state_cov = jacobian_t.matrix_transpose_multiply(tmp).as_scitbx_matrix()
      state_covs.append(state_cov)

    return state_covs

  def set_state_uncertainties(self, var_cov, multi_state_elt=None):
    """Send the calculated variance-covariance matrix for model state elements
    back to the model for storage alongside the model state, and potentially
    use in further propagation of error calculations.

    Args:
        var_cov: A variance-covariance matrix for the parameters of this model.
        multi_state_elt (bool): If the parameterisation is multi-state, a
            integer index is required to select which state this variance-
            covariance matrix refers to.
    """

    # To be implemented by the derived class, where it is clear what aspect
    # of the model under parameterisation is considered its state.
    pass
