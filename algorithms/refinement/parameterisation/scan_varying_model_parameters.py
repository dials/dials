#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

from __future__ import division
from dials.algorithms.refinement.parameterisation.model_parameters \
        import Parameter, ModelParameterisation
from math import exp
import abc
from scitbx import sparse

# reusable PHIL string for options affecting scan-varying parameterisation
phil_str = '''
smoother
  .help = "Options that affect scan-varying parameterisation"
  .expert_level = 1
{
  num_intervals = *fixed_width absolute
    .help = "Choose the way to determine the number of intervals for scan-"
            "varying refinement"
    .type = choice

  interval_width_degrees = 36.0
    .help = "Width of scan between checkpoints in degrees"
    .type = float(value_min=0.)

  absolute_num_intervals = 5
    .help = "Number of intervals between checkpoints if scan_varying"
            "refinement is requested"
    .type = int(value_min=1)
}
'''

class ScanVaryingParameterSet(Parameter):
  """Testing a class for a scan-varying parameter, in which values at rotation
  angle phi may be derived using smoothed interpolation between checkpoint
  values stored here. Externally, this is presented as a set of parameters.

  num_samples is the number of checkpoints. Other arguments are as Parameter.
  """

  def __init__(self, value, num_samples = 5, axis = None, ptype = None,
               name = "ScanVaryingParameterSet"):

    assert num_samples >= 2 #otherwise use scan-independent parameterisation

    value = [value] * num_samples
    self._name_stem = name
    name = [e + "_sample%d" % i for i, e in enumerate(
      [self._name_stem] * num_samples)]

    Parameter.__init__(self, value, axis, ptype, name)

    self._esd = [None] * num_samples
    self._num_samples = num_samples

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
  def name_stem(self):
    return self._name_stem

  def __str__(self):

    msg = "ScanVaryingParameterSet " + self.name_stem + ":\n"
    try:
      msg += "    Type: " + self.param_type + "\n"
    except TypeError:
      msg += "    Type: " + str(self.param_type) + "\n"
    try:
      msg += "    Axis: (%5.3f, %5.3f, %5.3f)" % tuple(self.axis) + "\n"
    except TypeError:
      msg += "    Axis: " + str(self.axis) + "\n"
    vals = ", ".join(["%5.3f"] * len(self)) % tuple(self.value)
    msg += "    Values: " + vals + "\n"
    try:
      sigs = ", ".join(["%5.3f"] * len(self)) % tuple(self.esd)
    except TypeError:
      sigs = ", ".join([str(e) for e in self.esd])
    msg += "    Sigmas: (" + sigs +") \n"

    return msg

class GaussianSmoother(object):
  """A Gaussian smoother for ScanVaryingModelParameterisations"""

  # Based largely on class SmoothedValue from Aimless.

  # Construct from range of raw unnormalised coordinate & number of sample
  # intervals. Set smoothing values to defaults, Nav = 3
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
      self._positions = [0.0, 1.0]
    elif self._nvalues == 3:
      self._positions = [0.0, 1.0, 2.0]
    else:
      self._positions = [e - 0.5 for e in range(self._nvalues)]

    # set default smoothing parameters
    self.set_smoothing(3, -1.0)

  def set_smoothing(self, num_average, sigma):
    """Set smoothing values:

    naverage: number of points included in each calculation
    sigma: width of the Gaussian used for smoothing.

    If sigma < 0, set to "optimum" (!) (or at least suitable) value from
    num_average """

    self._naverage = num_average
    if self._naverage > self._nvalues:
      self._naverage = self._nvalues
    self._half_naverage = self._naverage / 2.0
    self._sigma = sigma

    if self._naverage < 1 or self._naverage > 5:
      raise ValueError("num_average must be between 1 & 5")

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

    # use sparse storage as only self._naverage (default 3) values are non-zero
    weight = sparse.vector(self._nvalues)

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
        i1 = min(i1, self._nvalues - 2) # ensure separation of >= 2

    # now do stuff
    for i in range(i1, i2):

      ds = (z - self._positions[i]) / self._sigma
      w = exp(-ds*ds)
      sumwv += w * values[i]
      sumweight  += w
      weight[i] = w

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
  """Extending ModelParameterisation to deal with ScanVaryingParameterSets.

  For simplicity at this stage it is decreed that a
  ScanVaryingModelParameterisation consists only of ScanVaryingParameterSets.
  There is no combination with normal Parameters. This could be changed later,
  but there may be no reason to do so, hence starting with this simpler
  design"""

  # The initial state is here equivalent to the initial state of the
  # time static version of the parameterisation, as it is assumed that we
  # start with a flat model wrt rotation angle.

  __metaclass__  = abc.ABCMeta

  def __init__(self, model, initial_state, param_sets, smoother,
               experiment_ids, is_multi_state=False):

    ModelParameterisation.__init__(self, model, initial_state, param_sets,
      experiment_ids, is_multi_state)

    self._num_sets = len(self._param)
    self._set_len = len(param_sets[0])
    self._total_len = self._set_len * self._num_sets

    # ensure all internal parameter sets have the same number of parameters
    for param in self._param[1:]: assert len(param) == self._set_len

    # Link up with an object that will perform the smoothing.
    self._smoother = smoother
    assert self._smoother.num_values() == self._set_len

    # define an attribute for caching the variance-covariance matrix of
    # parameters
    self._var_cov = None

    return

  def num_free(self):
    """the number of free parameters"""

    if self._num_free is None:
      self._num_free = sum(not x.get_fixed() for x in self._param) * self._set_len
    return self._num_free

  # def num_total(self): inherited unchanged from ModelParameterisation

  def num_sets(self):
    """the number of parameter sets"""
    return self._num_sets

  @abc.abstractmethod
  def compose(self, t):
    """compose the model state at image number t from its initial state and
    its parameter list. Also calculate the derivatives of the state wrt
    each parameter in the list.

    Unlike ModelParameterisation, does not automatically update the actual
    model class. This should be done once refinement is complete."""

    pass

  def get_param_vals(self, only_free = True):
    """export the values of the internal list of parameters as a
    sequence of floats.

    If only_free, the values of fixed parameters are filtered from the
    returned list. Otherwise all parameter values are returned"""

    if only_free:
      return [x for e in self._param \
              if not e.get_fixed() for x in e.value]

    else:
      return [x for e in self._param for x in e.value]

  def get_param_names(self, only_free = True):
    """export the names of the internal list of parameters

    If only_free, the names of fixed parameters are filtered from the
    returned list. Otherwise all parameter names are returned"""

    # FIXME combine functionality with get_param_vals by returning a named,
    # ordered list?
    if only_free:
      return [x for e in self._param \
              if not e.get_fixed() for x in e.name]

    else:
      return [x for e in self._param for x in e.name]

  def set_param_vals(self, vals):
    """set the values of the internal list of parameters from a
    sequence of floats.

    First break the sequence into sub sequences of the same length
    as the _set_len.

    Only free parameter sets can have values assigned, therefore the
    length of vals must equal the value of num_free"""

    assert(len(vals) == self.num_free())
    i = 0
    for p in self._param:
      if not p.get_fixed(): # only set the free parameter sets
        new_vals = vals[i:i+self._set_len]
        p.value = new_vals
        i += self._set_len

    # compose with the new parameter values
    #self.compose()

    return

  def set_param_esds(self, esds):
    """set the estimated standard deviations of the internal list of parameters
    from a sequence of floats.

    First break the sequence into sub sequences of the same length
    as the _set_len.

    Only free parameters can be set, therefore the length of esds must equal
    the value of num_free"""

    assert(len(esds) == self.num_free())
    i = 0
    for p in self._param:
      if not p.get_fixed(): # only set the free parameter sets
        new_esds = esds[i:i+self._set_len]
        p.esd = new_esds
        i += self._set_len

    return

  #def get_fixed(self): inherited unchanged from ModelParameterisation

  #def set_fixed(self, fix): inherited unchanged from ModelParameterisation

  #def get_state(self): inherited unchanged from ModelParameterisation

  def get_ds_dp(self, only_free = True, use_none_as_null=False):
    """get a list of derivatives of the state wrt each parameter, as
    a list in the same order as the internal list of parameters. Requires
    compose to be called first at scan coordinate 't' so that each
    scan-dependent parameter is evaluated at coordinate t, corresponding to
    the original, unnormalised coordinates used to set up the smoother
    (t will most likely be along the dimension of image number).

    If only_free, the derivatives with respect to fixed parameters are
    omitted from the returned list. Otherwise a list for all parameters is
    returned, with null values for the fixed parameters.

    The internal list of derivatives self._dstate_dp may use None for null
    elements. By default these are converted to the null state, but
    optionally these may remain None to detect them easier and avoid
    doing calculations on null elements
    """

    if use_none_as_null:
      null = None
    else:
      null = self._null_state

    if only_free:
      return [null if ds_dp is None else ds_dp for row, p in zip(
              self._dstate_dp, self._param) if not p.get_fixed() for ds_dp in row]

    else:
      return [null if p.get_fixed() or ds_dp is None else ds_dp \
              for row, p in zip(self._dstate_dp, self._param) \
              for ds_dp in row]

  def get_smoothed_parameter_value(self, t, pset):
    """export the smoothed value of a parameter set at image number 't'
    using the smoother."""

    return self._smoother.value_weight(t, pset)[0]

  def calculate_state_uncertainties(self, var_cov=None):
    """Given a variance-covariance array for the parameters of this model,
    propagate those estimated errors into the uncertainties of the model state
    at every scan point"""

    if var_cov is not None:
      # first call, just cache the variance-covariance matrix
      self._var_cov = var_cov
      return None

    # later calls, make sure it has been cached! Otherwise ESDs cannot be
    # calculated, so return None
    if self._var_cov is None: return None

    # later calls, assumes compose has been called at image number t, so that
    # get_ds_dp will be specific for that image. Now call the base class method
    # and return the result
    return super(ScanVaryingModelParameterisation,
          self).calculate_state_uncertainties(self._var_cov)

  def set_state_uncertainties(self, var_cov_list):
    """Send the calculated variance-covariance matrices for model state elements
    for all scan points back to the model for storage alongside the model state
    """

    pass
