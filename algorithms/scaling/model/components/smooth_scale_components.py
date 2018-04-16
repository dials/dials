""""
Classes that each define a smoothly varying component of a scaling model.

These classes use a gaussian smoother (1D, 2D or 3D) to calculate the
inverse scale factors and derivatives with respect to the component
parameters.
"""
import numpy as np
from scitbx import sparse
from dials.array_family import flex
from dials_scaling_helpers_ext import row_multiply
from dials_refinement_helpers_ext import GaussianSmoother as GS1D
from dials_refinement_helpers_ext import GaussianSmoother2D as GS2D
from dials_refinement_helpers_ext import GaussianSmoother3D as GS3D
from dials_scratch_scaling_ext import elementwise_square
from dials.algorithms.scaling.model.components.scale_components import \
  ScaleComponentBase

# The following gaussian smoother classes make the implementation
# consistent with that used in dials.refinement.

class GaussianSmoother1D(GS1D):
  """A 1D Gaussian smoother."""

  def value_weight(self, x, value):
    """Return the value, weight and sumweight at a single point."""
    result = super(GaussianSmoother1D, self).value_weight(x,
      flex.double(value))
    return (result.get_value(), result.get_weight(), result.get_sumweight())

  def multi_value_weight(self, x, value):
    """Return the value, weight and sumweight at multiple points."""
    result = super(GaussianSmoother1D, self).multi_value_weight(
      flex.double(x),
      flex.double(value))
    return (result.get_value(), result.get_weight(), result.get_sumweight())

  def positions(self):
    """Return the smoother positions."""
    return list(super(GaussianSmoother1D, self).positions())

class GaussianSmoother2D(GS2D):
  """A 2D Gaussian smoother."""

  def value_weight(self, x, y, value):
    """Return the value, weight and sumweight at a single point."""
    result = super(GaussianSmoother2D, self).value_weight(x, y,
      flex.double(value))
    return (result.get_value(), result.get_weight(), result.get_sumweight())

  def multi_value_weight(self, x, y, value):
    """Return the value, weight and sumweight at multiple points."""
    result = super(GaussianSmoother2D, self).multi_value_weight(
      flex.double(x), flex.double(y),
      flex.double(value))
    return (result.get_value(), result.get_weight(), result.get_sumweight())

  def x_positions(self):
    """Return the smoother x-positions."""
    return list(super(GaussianSmoother2D, self).x_positions())

  def y_positions(self):
    """Return the smoother y-positions."""
    return list(super(GaussianSmoother2D, self).y_positions())

class GaussianSmoother3D(GS3D):
  """A 3D Gaussian smoother."""

  def value_weight(self, x, y, z, value):
    """Return the value, weight and sumweight at a single point."""
    result = super(GaussianSmoother3D, self).value_weight(x, y, z,
      flex.double(value))
    return (result.get_value(), result.get_weight(), result.get_sumweight())

  def multi_value_weight(self, x, y, z, value):
    """Return the value, weight and sumweight at multiple points."""
    result = super(GaussianSmoother3D, self).multi_value_weight(
      flex.double(x), flex.double(y), flex.double(z),
      flex.double(value))
    return (result.get_value(), result.get_weight(), result.get_sumweight())

  def x_positions(self):
    """Return the smoother x-positions."""
    return list(super(GaussianSmoother3D, self).x_positions())

  def y_positions(self):
    """Return the smoother y-positions."""
    return list(super(GaussianSmoother3D, self).y_positions())

  def z_positions(self):
    """Return the smoother z-positions."""
    return list(super(GaussianSmoother3D, self).z_positions())


class SmoothMixin(object):
  """Mixin class for smooth scale factor components.

  This uses a Gaussian smoother to calculate scales and derivatives
  based on the parameters and a have a set of normalised_values
  associated with the data."""

  def __init__(self):
    self._Vr = 1.0
    self._smoother = None

  @property
  def value(self):
    """Extra access to the parameters for the gaussian smoother."""
    return self._parameters

  @property
  def smoother(self):
    """The Gaussian smoother."""
    return self._smoother

  @staticmethod
  def nparam_to_val(n_params):
    """Convert the number of parameters to the required input value
    for the smoother."""
    assert n_params >= 2, '''cannot initialise a smooth scale factor
      for <2 parameters.'''
    if n_params == 2 or n_params == 3:
      return n_params - 1
    return n_params - 2

class SmoothScaleComponent1D(ScaleComponentBase, SmoothMixin):
  """A smoothly varying scale component in one dimension."""

  def __init__(self, initial_values, col_name, parameter_esds=None):
    super(SmoothScaleComponent1D, self).__init__(initial_values,
      parameter_esds)
    self._normalised_values = None
    self._col_name = col_name

  @property
  def normalised_values(self):
    """This is a list of the relevant data needed to calculate the
    inverse scale factors, normalised to give 'normalised coordinate
    values' that fit in the range of the smoother parameters, which
    are defined as a 1D array at normalised coordinates separated by
    a spacing of 1."""
    return self._normalised_values

  @property
  def col_name(self):
    """The column name to use to obtain normalised coordinates from a
    reflection table."""
    return self._col_name

  def update_reflection_data(self, reflection_table, selection=None):
    """Set the normalised coordinate values and configure the smoother."""
    normalised_values = reflection_table[self._col_name]
    if selection:
      normalised_values = normalised_values.select(selection)
    # Make sure zeroed correctly.
    normalised_values = normalised_values - min(normalised_values)
    self._normalised_values = normalised_values
    phi_range_deg = [int(min(self._normalised_values)//1),
                     int(max(self._normalised_values)//1)+1]
    self._smoother = GaussianSmoother1D(phi_range_deg,
      self.nparam_to_val(self._n_params))
    self.inverse_scales = flex.double(normalised_values.size(), 1.0)
    self._n_refl = self.inverse_scales.size()

  def calculate_scales_and_derivatives(self, curvatures=False):
    value, weight, sumweight = self._smoother.multi_value_weight(
      self._normalised_values, self.value)
    inv_sw = 1. / sumweight
    dv_dp = row_multiply(weight, inv_sw)
    self._inverse_scales = value
    self._derivatives = dv_dp
    if curvatures:
      self._curvatures = sparse.matrix(self._inverse_scales.size(), self.n_params)

  def calculate_scales(self):
    """"Only calculate the scales if needed, for performance."""
    value, _, _ = self._smoother.multi_value_weight(self._normalised_values,
      self.value)
    self._inverse_scales = value


class SmoothBScaleComponent1D(SmoothScaleComponent1D):
  '''Subclass of SmoothScaleComponent1D to implement a smoothly
  varying B-factor correction.'''

  def __init__(self, initial_values, col_name, parameter_esds=None):
    super(SmoothBScaleComponent1D, self).__init__(initial_values, col_name,
      parameter_esds)
    self._d_values = None

  @property
  def d_values(self):
    """The current set of d-values associated with this component."""
    return self._d_values

  def update_reflection_data(self, reflection_table, selection=None):
    super(SmoothBScaleComponent1D, self).update_reflection_data(
      reflection_table, selection)
    self._d_values = reflection_table['d']
    if selection:
      self._d_values = self._d_values.select(selection)

  def calculate_scales_and_derivatives(self, curvatures=False):
    super(SmoothBScaleComponent1D, self).calculate_scales_and_derivatives()
    self._inverse_scales = flex.double(np.exp(self._inverse_scales
      /(2.0 * (self._d_values**2))))
    self._derivatives = row_multiply(self._derivatives,
      self._inverse_scales / (2.0 * (self._d_values**2)))
    if curvatures:
      self._curvatures = row_multiply(elementwise_square(self._derivatives),
        1.0/self._inverse_scales)

  def calculate_scales(self):
    super(SmoothBScaleComponent1D, self).calculate_scales()
    self._inverse_scales = flex.double(np.exp(self._inverse_scales
      /(2.0 * (self._d_values**2))))


class SmoothScaleComponent2D(ScaleComponentBase, SmoothMixin):
  """Implementation of a 2D array-based smoothly varying scale factor.

  A 2d array of parameters is defined, and the scale factor at fractional
  coordinates is calculated as smoothly varying based on the distance to
  the nearby parameters as calculated in the GaussianSmoother2D. The
  initial values are passed as a 1D array, and shape is a 2-tuple
  indicating the number of parameters in each dimension."""

  def __init__(self, initial_values, shape, col_names, parameter_esds=None):
    assert len(initial_values) == (shape[0] * shape[1]), '''The shape
    information to initialise a 2D smoother is inconsistent with the length
    of the initial parameter list.'''
    super(SmoothScaleComponent2D, self).__init__(initial_values, parameter_esds)
    self._n_x_params = shape[0]
    self._n_y_params = shape[1]
    self._col_names = col_names
    self._normalised_x_values = None
    self._normalised_y_values = None

  @property
  def col_names(self):
    """The column names used to obtain normalised coordinates from a
    reflection table."""
    return self._col_names

  @property
  def n_x_params(self):
    """The number of parameters that parameterise the x-component."""
    return self._n_x_params

  @property
  def n_y_params(self):
    """The number of parameters that parameterise the y-component."""
    return self._n_y_params

  @property
  def normalised_x_values(self):
    """The normalised coordinate values in the first dimension."""
    return self._normalised_x_values

  @property
  def normalised_y_values(self):
    """The normalised coordinate values in the second dimension."""
    return self._normalised_y_values

  def update_reflection_data(self, reflection_table, selection=None):
    '''control access to setting all of reflection data at once'''
    if selection:
      reflection_table = reflection_table.select(selection)
    normalised_x_values = reflection_table[self._col_names[0]]
    normalised_y_values = reflection_table[self._col_names[1]]
    self._normalised_x_values = normalised_x_values - min(normalised_x_values)
    self._normalised_y_values = normalised_y_values - min(normalised_y_values)
    x_range = [int(min(self._normalised_x_values)//1),
               int(max(self._normalised_x_values)//1)+1]
    y_range = [int(min(self._normalised_y_values)//1),
               int(max(self._normalised_y_values)//1)+1]
    self._smoother = GaussianSmoother2D(x_range, self.nparam_to_val(
      self._n_x_params), y_range, self.nparam_to_val(self._n_y_params))
    self.inverse_scales = flex.double(self._normalised_x_values.size(), 1.0)
    self._n_refl = self.inverse_scales.size()

  def calculate_scales_and_derivatives(self, curvatures=False):
    value, weight, sumweight = self._smoother.multi_value_weight(
      self._normalised_x_values, self._normalised_y_values, self.value)
    inv_sw = 1. / sumweight
    dv_dp = row_multiply(weight, inv_sw)
    self._inverse_scales = value
    self._derivatives = dv_dp
    if curvatures:
      self._curvatures = sparse.matrix(self._inverse_scales.size(), self.n_params)

  def calculate_scales(self):
    """"Only calculate the scales if needed, for performance."""
    value, _, _ = self._smoother.multi_value_weight(self._normalised_x_values,
      self._normalised_y_values, self.value)
    self._inverse_scales = value


class SmoothScaleComponent3D(ScaleComponentBase, SmoothMixin):
  """Implementation of a 3D array-based smoothly varying scale factor.

  A 3d array of parameters is defined, and the scale factor at fractional
  coordinates is calculated as smoothly varying based on the distance to
  the nearby parameters as calculated in the GaussianSmoother3D. The
  initial values are passed as a 1D array, and shape is a 3-tuple
  indicating the number of parameters in each dimension."""

  def __init__(self, initial_values, shape, col_names, parameter_esds=None):
    assert len(initial_values) == (shape[0] * shape[1] * shape[2]), '''The
    shape information to initialise a 3D smoother is inconsistent with the
    length of the initial parameter list.'''
    super(SmoothScaleComponent3D, self).__init__(initial_values,
      parameter_esds)
    self._n_x_params = shape[0]
    self._n_y_params = shape[1]
    self._n_z_params = shape[2]
    self._normalised_x_values = None
    self._normalised_y_values = None
    self._normalised_z_values = None
    self._col_names = col_names

  @property
  def col_names(self):
    """The column names used to obtain normalised coordinates from a
    reflection table."""
    return self._col_names

  @property
  def n_x_params(self):
    """The number of parameters that parameterise the x-component."""
    return self._n_x_params

  @property
  def n_y_params(self):
    """The number of parameters that parameterise the y-component."""
    return self._n_y_params

  @property
  def n_z_params(self):
    """The number of parameters that parameterise the z-component."""
    return self._n_z_params

  @property
  def normalised_x_values(self):
    """The normalised coordinate values in the first dimension."""
    return self._normalised_x_values

  @property
  def normalised_y_values(self):
    """The normalised coordinate values in the second dimension."""
    return self._normalised_y_values

  @property
  def normalised_z_values(self):
    """The normalised coordinate values in the third dimension."""
    return self._normalised_z_values

  def update_reflection_data(self, reflection_table, selection=None):
    '''control access to setting all of reflection data at once'''
    if selection:
      reflection_table = reflection_table.select(selection)
    normalised_x_values = reflection_table[self._col_names[0]]
    normalised_y_values = reflection_table[self._col_names[1]]
    normalised_z_values = reflection_table[self._col_names[2]]
    """Set the normalised coordinate values and configure the smoother."""
    self._normalised_x_values = normalised_x_values - min(normalised_x_values)
    self._normalised_y_values = normalised_y_values - min(normalised_y_values)
    self._normalised_z_values = normalised_z_values - min(normalised_z_values)
    x_range = [int(min(self._normalised_x_values)//1),
               int(max(self._normalised_x_values)//1)+1]
    y_range = [int(min(self._normalised_y_values)//1),
               int(max(self._normalised_y_values)//1)+1]
    z_range = [int(min(self._normalised_z_values)//1),
               int(max(self._normalised_z_values)//1)+1]
    self._smoother = GaussianSmoother3D(x_range, self.nparam_to_val(
      self._n_x_params), y_range, self.nparam_to_val(self._n_y_params),
      z_range, self.nparam_to_val(self._n_z_params))
    self.inverse_scales = flex.double(self._normalised_x_values.size(), 1.0)
    self._n_refl = self.inverse_scales.size()

  def calculate_scales_and_derivatives(self, curvatures=False):
    value, weight, sumweight = self._smoother.multi_value_weight(
      self._normalised_x_values, self._normalised_y_values,
      self._normalised_z_values, self.value)
    inv_sw = 1. / sumweight
    dv_dp = row_multiply(weight, inv_sw)
    self._inverse_scales = value
    self._derivatives = dv_dp
    if curvatures:
      self._curvatures = sparse.matrix(self._inverse_scales.size(), self.n_params)

  def calculate_scales(self):
    """"Only calculate the scales if needed, for performance."""
    value, _, _ = self._smoother.multi_value_weight(self._normalised_x_values,
      self._normalised_y_values, self._normalised_z_values, self.value)
    self._inverse_scales = value
