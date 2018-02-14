'''
These classes define objects to hold a flex array of scale factors,
including the required methods for calculating smooth scale factors and
'''
import abc
import numpy as np
from dials.array_family import flex
from dials.algorithms.refinement.parameterisation.scan_varying_model_parameters \
        import GaussianSmoother#, GaussianSmoother2D, GaussianSmoother3D
from dials_scaling_helpers_ext import row_multiply
from scitbx import sparse
from dials_refinement_helpers_ext import GaussianSmoother2D as GS2D
from dials_refinement_helpers_ext import GaussianSmoother3D as GS3D
#from dials_refinement_helpers_ext import GaussianSmoother2D as GS2D
#from dials.algorithms.refinement.boost_python.gaussian_smoother_2D import \

class GaussianSmoother2D(GS2D):
  """A 2D Gaussian smoother"""

  def value_weight(self, x, y, param):
    result = super(GaussianSmoother2D, self).value_weight(x, y,
      flex.double(param.value))
    return (result.get_value(), result.get_weight(), result.get_sumweight())

  def multi_value_weight(self, x, y, param):
    result = super(GaussianSmoother2D, self).multi_value_weight(
      flex.double(x), flex.double(y),
      flex.double(param.value))
    return (result.get_value(), result.get_weight(), result.get_sumweight())

  def x_positions(self):
    return list(super(GaussianSmoother2D, self).x_positions())

  def y_positions(self):
    return list(super(GaussianSmoother2D, self).y_positions())

class GaussianSmoother3D(GS3D):
  """A 3D Gaussian smoother"""

  def value_weight(self, x, y, z, param):
    result = super(GaussianSmoother3D, self).value_weight(x, y, z,
      flex.double(param.value))
    return (result.get_value(), result.get_weight(), result.get_sumweight())

  def multi_value_weight(self, x, y, z, param):
    result = super(GaussianSmoother3D, self).multi_value_weight(
      flex.double(x), flex.double(y), flex.double(z),
      flex.double(param.value))
    return (result.get_value(), result.get_weight(), result.get_sumweight())

  def x_positions(self):
    return list(super(GaussianSmoother3D, self).x_positions())

  def y_positions(self):
    return list(super(GaussianSmoother3D, self).y_positions())

  def z_positions(self):
    return list(super(GaussianSmoother3D, self).z_positions())




class ScaleFactor(object):
  '''Base ScaleFactor class, with an interface to access the parameters,
  inverse_scales and derivatives.'''
  __metaclass__ = abc.ABCMeta

  def __init__(self, initial_value, scaling_options=None):
    #if isinstance(initial_value, list):
    #  assert len(initial_value) == n_parameters
    #  self._parameters = flex.double(initial_value)
    #else:
    self._parameters = initial_value#flex.double([initial_value] * n_parameters)
    self._n_params = len(self._parameters)
    self._scaling_options = scaling_options
    self._inverse_scales = None
    self._derivatives = None

  @property
  def n_params(self):
    '''number of parameters that parameterise the scale factor'''
    return self._n_params

  @property
  def parameters(self):
    '''parameters of the model'''
    return self._parameters

  @parameters.setter
  def parameters(self, values):
    if len(values) != len(self._parameters):
      assert 0, '''attempting to set a new set of parameters of different
      length than previous assignment: was %s, attempting %s''' % (
        len(self._parameters), len(values))
    self._parameters = values

  @property
  def inverse_scales(self):
    '''inverse scale factors associated with reflections'''
    return self._inverse_scales

  @inverse_scales.setter
  def inverse_scales(self, new_inverse_scales):
    self._inverse_scales = new_inverse_scales

  @property
  def derivatives(self):
    '''derivatives matrix of reflections with respect to model parameters'''
    return self._derivatives

  @abc.abstractmethod
  def update_reflection_data(self):
    '''method to be filled in by subclasses. Force setting of normalised
    values and any other data (e.g. d_values) at same time.'''
    pass

  def calculate_scales_and_derivatives(self):
    """method to be filled in by subclasses"""
    pass


class KScaleFactor(ScaleFactor):
  '''ScaleFactor object for a single global scale parameter.'''
  def __init__(self, initial_value, scaling_options=None):
    assert len(initial_value) == 1, '''Can only have one parameter for a K or B
      Scale Factor'''
    super(KScaleFactor, self).__init__(initial_value, scaling_options)
    self._n_refl = None

  @property
  def n_refl(self):
    return self._n_refl

  def update_reflection_data(self, n_refl):
    self._n_refl = n_refl

  def calculate_scales_and_derivatives(self):
    self._inverse_scales = flex.double([self._parameters[0]] * self.n_refl)
    self._derivatives = sparse.matrix(self.n_refl, 1)
    for i in range(self.n_refl):
      self._derivatives[i, 0] = 1.0


class BScaleFactor(KScaleFactor):
  '''ScaleFactor object for a single global B-scale parameter.'''
  def __init__(self, initial_value, scaling_options=None):
    super(BScaleFactor, self).__init__(initial_value, scaling_options)
    self._d_values = None

  @property
  def d_values(self):
    return self._d_values

  def update_reflection_data(self, dvalues):
    '''also sets n_refl. Allow to set to a different length to previous.'''
    self._d_values = dvalues
    self._n_refl = len(dvalues)

  def calculate_scales_and_derivatives(self):
    self._inverse_scales = flex.double(np.exp(flex.double(
      [self._parameters[0]] * self._n_refl) / (2.0 * (self._d_values**2))))
    self._derivatives = sparse.matrix(self._n_refl, 1)
    for i in range(self._n_refl):
      self._derivatives[i, 0] = (self._inverse_scales[i]
        / (2.0 * (self._d_values[i]**2)))

class SHScaleFactor(ScaleFactor):
  '''ScaleFactor class for a spherical harmonic absorption correction'''
  def __init__(self, initial_value, scaling_options=None):
    super(SHScaleFactor, self).__init__(initial_value, scaling_options)
    self._harmonic_values = None

  @property
  def harmonic_values(self):
    return self._harmonic_values

  def update_reflection_data(self, harmonic_values):
    self._harmonic_values = harmonic_values
    self.calculate_scales_and_derivatives()

  def calculate_scales_and_derivatives(self):
    '''calculation of scale factors and derivatives from
       spherical harmonic coefficient matrix'''
    abs_scale = flex.double([1.0] * self._harmonic_values.n_rows)#unity term
    for i, col in enumerate(self._harmonic_values.cols()):
      abs_scale += flex.double(col.as_dense_vector() * self._parameters[i])
    self._inverse_scales = abs_scale
    self._derivatives = self._harmonic_values

class SmoothScaleFactor(ScaleFactor):
  '''Base class for Smooth ScaleFactor objects - which allow use of a
  Gaussian smoother to calculate scales and derivatives based on the
  parameters and a set of normalised_values associated with each datapoint.'''
  def __init__(self, initial_value, scaling_options=None):
    super(SmoothScaleFactor, self).__init__(initial_value, scaling_options)
    self._normalised_values = None
    self.smoothing_window = 2.5 #redundant - replace with link to num_average
    self.Vr = 1.0 #replace with link to sigma of gaussian smoother?
    self._smoother = None #placeholder for gaussian smoother
    self.normalisation_interval = None

  @property
  def value(self):
    '''add extra access for gaussian smoother'''
    return self._parameters

  @property
  def normalised_values(self):
    '''normalised_values is the column from the reflection table
    of the normalised time/resolution etc'''
    return self._normalised_values

  def update_reflection_data(self):
    pass

class SmoothScaleFactor1D(SmoothScaleFactor):
  '''Class to implement a smooth scale factor in one dimension.'''

  def update_reflection_data(self, normalised_values):
    '''control access to setting all of reflection data at once'''
    self._normalised_values = normalised_values
    phi_range_deg = [int(min(self._normalised_values)//1),
                     int(max(self._normalised_values)//1)+1]
    #include an assertion here to stop setting too few params?
    self._smoother = GaussianSmoother(phi_range_deg, self._n_params - 2)
    self.inverse_scales = flex.double([1.0]*len(normalised_values))

  def calculate_scales_and_derivatives(self):
    value, weight, sumweight = self._smoother.multi_value_weight(
      self._normalised_values, self)
    inv_sw = 1. / sumweight
    dv_dp = row_multiply(weight, inv_sw)
    self._inverse_scales = value
    self._derivatives = dv_dp

  def calculate_scales(self):
    '''method to only calculate scales if this is needed, for performance.'''
    value, _, _ = self._smoother.multi_value_weight(self._normalised_values, self)
    self._inverse_scales = value

class SmoothBScaleFactor1D(SmoothScaleFactor1D):
  '''Subclass to SmoothScaleFactor1D for a smooth B-scale correction.'''
  def __init__(self, initial_value, scaling_options=None):
    super(SmoothBScaleFactor1D, self).__init__(initial_value, scaling_options)
    #self.Vr = 0.5
    self._d_values = None

  @property
  def d_values(self):
    return self._d_values

  def update_reflection_data(self, normalised_values, dvalues):
    assert len(normalised_values) == len(dvalues)
    super(SmoothBScaleFactor1D, self).update_reflection_data(normalised_values)
    self._d_values = dvalues

  def calculate_scales_and_derivatives(self):
    super(SmoothBScaleFactor1D, self).calculate_scales_and_derivatives()
    self._inverse_scales = flex.double(np.exp(self._inverse_scales
      /(2.0 * (self._d_values**2))))
    self._derivatives = row_multiply(self._derivatives,
      self._inverse_scales / (2.0 * (self._d_values**2)))

  def calculate_scales(self):
    super(SmoothBScaleFactor1D, self).calculate_scales()
    self._inverse_scales = flex.double(np.exp(self._inverse_scales
      /(2.0 * (self._d_values**2))))




class SmoothScaleFactor2D(SmoothScaleFactor):
  '''Class to implement a smooth scale factor in two dimensions.'''
  def __init__(self, initial_value, shape, scaling_options=None):
    super(SmoothScaleFactor2D, self).__init__(initial_value, scaling_options)
    self.n_x_params = shape[0]
    self.n_y_params = shape[1]
    self._normalised_x_values = None
    self._normalised_y_values = None

  @property
  def normalised_x_values(self):
    return self._normalised_x_values

  @property
  def normalised_y_values(self):
    return self._normalised_y_values

  def update_reflection_data(self, normalised_x_values, normalised_y_values):
    '''control access to setting all of reflection data at once'''
    self._normalised_x_values = normalised_x_values
    self._normalised_y_values = normalised_y_values
    x_range = [int(min(self._normalised_x_values)//1),
               int(max(self._normalised_x_values)//1)+1]
    y_range = [int(min(self._normalised_y_values)//1),
               int(max(self._normalised_y_values)//1)+1]
    #include an assertion here to stop setting too few params?
    self._smoother = GaussianSmoother2D(x_range, self.n_x_params - 2,
      y_range, self.n_y_params - 2)
    self.inverse_scales = flex.double([1.0]*len(normalised_x_values))

  def calculate_scales_and_derivatives(self):
    value, weight, sumweight = self._smoother.multi_value_weight(
      self._normalised_x_values, self._normalised_y_values, self)
    inv_sw = 1. / sumweight
    dv_dp = row_multiply(weight, inv_sw)
    self._inverse_scales = value
    self._derivatives = dv_dp

  def calculate_scales(self):
    '''method to only calculate scales if this is needed, for performance.'''
    value, _, _ = self._smoother.multi_value_weight(self._normalised_x_values,
      self._normalised_y_values, self)
    self._inverse_scales = value

class SmoothScaleFactor3D(SmoothScaleFactor2D):
  '''Class to implement a smooth scale factor in two dimensions.'''
  def __init__(self, initial_value, shape, scaling_options=None):
    super(SmoothScaleFactor3D, self).__init__(initial_value, shape, scaling_options)
    self.n_z_params = shape[2]
    self._normalised_z_values = None

  @property
  def normalised_z_values(self):
    return self._normalised_z_values

  def update_reflection_data(self, normalised_x_values, normalised_y_values,
      normalised_z_values):
    '''control access to setting all of reflection data at once'''
    self._normalised_x_values = normalised_x_values
    self._normalised_y_values = normalised_y_values
    self._normalised_z_values = normalised_z_values
    x_range = [int(min(self._normalised_x_values)//1),
               int(max(self._normalised_x_values)//1)+1]
    y_range = [int(min(self._normalised_y_values)//1),
               int(max(self._normalised_y_values)//1)+1]
    z_range = [int(min(self._normalised_z_values)//1),
               int(max(self._normalised_z_values)//1)+1]
    #include an assertion here to stop setting too few params?
    self._smoother = GaussianSmoother3D(x_range, self.n_x_params - 2,
      y_range, self.n_y_params - 2, z_range, self.n_z_params - 2)
    self.inverse_scales = flex.double([1.0]*len(normalised_x_values))

  def calculate_scales_and_derivatives(self):
    value, weight, sumweight = self._smoother.multi_value_weight(
      self._normalised_x_values, self._normalised_y_values,
      self._normalised_z_values, self)
    inv_sw = 1. / sumweight
    dv_dp = row_multiply(weight, inv_sw)
    self._inverse_scales = value
    self._derivatives = dv_dp

  def calculate_scales(self):
    '''method to only calculate scales if this is needed, for performance.'''
    value, _, _ = self._smoother.multi_value_weight(self._normalised_x_values,
      self._normalised_y_values, self._normalised_z_values, self)
    self._inverse_scales = value
