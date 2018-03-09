""""
Classes that each define one component of a scaling model.

These classes hold the component parameters and a component of
the inverse scale factors, as well as methods to calculate and
set the inverse scale factors and derivatives.
"""
import abc
import numpy as np
from dials.array_family import flex
from scitbx import sparse


class ScaleComponentBase(object):
  """Base scale component class.

  This defines an interface to access the parameters, the component
  of the inverse scale factor and it's derivatives with respect to
  the parameters. Scale components derived from the base class are
  designed to be instantiated by a ScalingModel class, by supplying
  an initial array of parameters and optionally the current estimated
  standard deviations. The relevant data from a reflection table is
  added later by a Scaler using the update_reflection_data method.
  This behaviour allows data to easily be added/changed after selecting
  subsets of the data."""

  __metaclass__ = abc.ABCMeta

  def __init__(self, initial_values, parameter_esds=None):
    self._parameters = initial_values
    self._parameter_esds = parameter_esds
    self._var_cov = None
    self._n_params = len(self._parameters)
    self._inverse_scales = None
    self._inverse_scale_variances = None
    self._derivatives = None
    self._curvatures = 0.0

  @property
  def n_params(self):
    """The number of parameters that parameterise the component."""
    return self._n_params

  @property
  def parameters(self):
    """The parameters of the component."""
    return self._parameters

  @parameters.setter
  def parameters(self, values):
    if len(values) != len(self._parameters):
      assert 0, '''attempting to set a new set of parameters of different
      length than previous assignment: was %s, attempting %s''' % (
        len(self._parameters), len(values))
    self._parameters = values

  @property
  def parameter_esds(self):
    """The estimated standard deviations of the parameters."""
    return self._parameter_esds

  @parameter_esds.setter
  def parameter_esds(self, esds):
    assert len(esds) == len(self._parameters)
    self._parameter_esds = esds

  @property
  def var_cov_matrix(self):
    """The variance/covariance matrix of the parameters."""
    return self._var_cov

  @var_cov_matrix.setter
  def var_cov_matrix(self, var_cov):
    self._var_cov = var_cov

  @property
  def inverse_scales(self):
    """The inverse scale factors associated with the data."""
    return self._inverse_scales

  @inverse_scales.setter
  def inverse_scales(self, new_inverse_scales):
    self._inverse_scales = new_inverse_scales

  @property
  def derivatives(self):
    """A spares matrix of the derivatives of the inverse scale
    factors with respect to the component parameters."""
    return self._derivatives

  @property
  def curvatures(self):
    """A spares matrix of the curvatures of the inverse scale
    factors with respect to the component parameters."""
    return self._curvatures

  @abc.abstractmethod
  def update_reflection_data(self):
    """Add or change the relevant reflection data for the component.

    No restrictions should be placed on the data remaining the same
    size, to allow easy changing of the data contained during scaling.
    The input data will be specific to the component."""
    pass

  @abc.abstractmethod
  def calculate_scales_and_derivatives(self):
    """Use the component parameters to calculate and set
    self._inverse_scales and self._derivatives."""
    pass


class SingleScaleFactor(ScaleComponentBase):
  """A model component consisting of a single global scale parameter.

  The inverse scale factor for every reflection is the parameter
  value and the derivatives are all 1.0."""

  def __init__(self, initial_values, parameter_esds=None):
    assert len(initial_values) == 1, '''This model component is only designed
      for a single global scale component'''
    super(SingleScaleFactor, self).__init__(initial_values, parameter_esds)
    self._n_refl = None

  @property
  def n_refl(self):
    """The current number of reflections associated with this component."""
    return self._n_refl

  def update_reflection_data(self, n_refl):
    """Add reflection data to the component, only n_reflections needed."""
    self._n_refl = n_refl

  def calculate_scales_and_derivatives(self):
    self._inverse_scales = flex.double([self._parameters[0]] * self.n_refl)
    self._derivatives = sparse.matrix(self.n_refl, 1)
    for i in range(self.n_refl):
      self._derivatives[i, 0] = 1.0


class SingleBScaleFactor(ScaleComponentBase):
  """A model component for a single global B-factor parameter.

  The inverse scale factor for each reflection is given by
  S = exp(B/(2 * d^2)), the derivatives are S/(2 * d^2)."""

  def __init__(self, initial_values, parameter_esds=None):
    super(SingleBScaleFactor, self).__init__(initial_values, parameter_esds)
    self._d_values = None
    self._n_refl = None

  @property
  def d_values(self):
    """The current set of d-values associated with this component."""
    return self._d_values

  @property
  def n_refl(self):
    """The current number of reflections associated with this component."""
    return self._n_refl

  def update_reflection_data(self, dvalues):
    """"Add reflection data to the component, only the d-values and number
    of reflections are needed."""
    self._d_values = dvalues
    self._n_refl = len(dvalues)

  def calculate_scales_and_derivatives(self):
    self._inverse_scales = flex.double(np.exp(flex.double(
      [self._parameters[0]] * self._n_refl) / (2.0 * (self._d_values**2))))
    self._derivatives = sparse.matrix(self._n_refl, 1)
    for i in range(self._n_refl):
      self._derivatives[i, 0] = (self._inverse_scales[i]
        / (2.0 * (self._d_values[i]**2)))


class SHScaleComponent(ScaleComponentBase):
  """A model component for a spherical harmonic absorption correction.

  This component uses a set of spherical harmonic functions to define
  an absorption surface for the crystal. A matrix of spherical harmonic
  coefficients for the data is stored in self._harmonic_values and is
  used to calculate the scales and derivatives.
  The scale is given by S = 1 + (sum_l sum_m Clm * Ylm) where Clm are
  the model parameters and Ylm are the spherical harmonic coefficients,
  the derivatives are then simply the coefficients Ylm."""

  def __init__(self, initial_values, parameter_esds=None):
    super(SHScaleComponent, self).__init__(initial_values, parameter_esds)
    self._harmonic_values = None

  @property
  def harmonic_values(self):
    """A matrix of harmonic coefficients for the data."""
    return self._harmonic_values

  def update_reflection_data(self, harmonic_values):
    """Update the spherical harmonic coefficients."""
    self._harmonic_values = harmonic_values
    self.calculate_scales_and_derivatives()

  def calculate_scales_and_derivatives(self):
    abs_scale = flex.double([1.0] * self._harmonic_values.n_rows) # Unity term
    for i, col in enumerate(self._harmonic_values.cols()):
      abs_scale += flex.double(col.as_dense_vector() * self._parameters[i])
    self._inverse_scales = abs_scale
    self._derivatives = self._harmonic_values
