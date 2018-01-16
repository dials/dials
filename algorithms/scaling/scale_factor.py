'''
These classes define objects to hold a flex array of scale factors,
including the required methods for calculating smooth scale factors and
'''
import abc
from dials.array_family import flex
import numpy as np
#import math as math
from dials.algorithms.refinement.parameterisation.scan_varying_model_parameters \
        import GaussianSmoother
from dials_scaling_helpers_ext import row_multiply
from scitbx import sparse

class ScaleFactor(object):
  '''Base ScaleFactor class, with an interface to access the parameters,
  inverse_scales and derivatives.'''
  __metaclass__ = abc.ABCMeta

  def __init__(self, initial_value, n_parameters, scaling_options=None):
    self._parameters = flex.double([initial_value] * n_parameters)
    self._n_params = n_parameters
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
    super(KScaleFactor, self).__init__(initial_value, 1, scaling_options)
    self._n_refl = None#n_refl

  @property
  def n_refl(self):
    return self._n_refl

  def update_reflection_data(self, n_refl):
    self._n_refl = n_refl

  #@n_refl.setter
  #def n_refl(self, number_of_refl):
  #  self._n_refl = number_of_refl

  def calculate_scales_and_derivatives(self):
    self._inverse_scales = flex.double([self._parameters[0]] * self.n_refl)
    self._derivatives = sparse.matrix(self.n_refl, 1)
    for i in range(self.n_refl):
      self._derivatives[i, 0] = 1.0


class BScaleFactor(ScaleFactor):
  '''ScaleFactor object for a single global B-scale parameter.'''
  def __init__(self, initial_value, scaling_options=None):
    super(BScaleFactor, self).__init__(initial_value, 1, scaling_options)
    self._d_values = None#d_values
    self._n_refl = None#len(d_values)

  @property
  def d_values(self):
    return self._d_values

  def update_reflection_data(self, dvalues):
    '''also sets n_refl. Allow to set to a different length to previous.'''
    self._d_values = dvalues
    self._n_refl = len(dvalues)
    #self._n_refl = number_of_refl

  #@d_values.setter
  #def d_values(self, new_values):
  #  '''also sets n_refl. Allow to set to a different length to previous.'''
  #  self._d_values = new_values
  #  self._n_refl = len(new_values)

  def calculate_scales_and_derivatives(self):
    self._inverse_scales = flex.double(np.exp(flex.double(
      [self._parameters[0]] * self._n_refl) / (2.0 * (self._d_values**2))))
    self._derivatives = sparse.matrix(self._n_refl, 1)
    for i in range(self._n_refl):
      self._derivatives[i, 0] = (self._inverse_scales[i]
                                 / (2.0 * (self._d_values[i]**2)))


class SmoothScaleFactor(ScaleFactor):
  '''Base class for Smooth ScaleFactor objects - which allow use of a
  Gaussian smoother to calculate scales and derivatives based on the
  parameters and a set of normalised_values associated with each datapoint.'''
  def __init__(self, initial_value, n_parameters, scaling_options=None):
    super(SmoothScaleFactor, self).__init__(initial_value, n_parameters,
      scaling_options)
    self._normalised_values = None
    self.Vr = 1.0# 0.7

  @property
  def value(self):
    '''add extra access for gaussian smoother'''
    return self._parameters

  @property
  def normalised_values(self):
    '''normalised_values is the column from the reflection table
    of the normalised time/resolution etc'''
    return self._normalised_values


class SmoothScaleFactor1D(SmoothScaleFactor):
  '''Class to implement a smooth scale factor in one dimension.'''
  def __init__(self, initial_value, n_parameters, scaling_options=None):
    super(SmoothScaleFactor1D, self).__init__(initial_value, n_parameters,
      scaling_options)
    self.smoothing_window = 2.5 #must be less than 3 to avoid indexing errors
    self._smoother = None #placeholder for gaussian smoother

  def update_reflection_data(self, normalised_values):
    '''control access to setting all of reflection data at once'''
    self._normalised_values = normalised_values
    phi_range_deg = [int(min(self._normalised_values)//1),
                     int(max(self._normalised_values)//1)+1]
    self._smoother = GaussianSmoother(phi_range_deg, self._n_params - 2)
    #self.inverse_scales = flex.double([1.0]*len(new_values))
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
  def __init__(self, initial_value, n_parameters, scaling_options=None):
    super(SmoothBScaleFactor1D, self).__init__(initial_value, n_parameters,
      scaling_options)
    self.Vr = 0.5
    self._d_values = None
    self._n_refl = None

  @property
  def d_values(self):
    return self._d_values

  def update_reflection_data(self, normalised_values, dvalues):
    assert len(normalised_values) == len(dvalues)
    self._normalised_values = normalised_values
    phi_range_deg = [int(min(self._normalised_values)//1),
                     int(max(self._normalised_values)//1)+1]
    self._smoother = GaussianSmoother(phi_range_deg, self._n_params - 2)
    self.inverse_scales = flex.double([1.0]*len(normalised_values))
    self._d_values = dvalues
    self._n_refl = len(dvalues)

  def calculate_scales_and_derivatives(self):
    value, weight, sumweight = self._smoother.multi_value_weight(
      self._normalised_values, self)
    self._inverse_scales = flex.double(np.exp(value/(2.0 * (self._d_values**2))))
    inv_sw = 1. / sumweight
    dv_dp = row_multiply(weight, inv_sw)
    multiplicative_factor = self._inverse_scales / (2.0 * (self._d_values**2))
    self._derivatives = row_multiply(dv_dp, multiplicative_factor)

  def calculate_scales(self):
    value, _, _ = self._smoother.multi_value_weight(self._normalised_values, self)
    self._inverse_scales = flex.double(np.exp(value/(2.0 * (self._d_values**2))))


class SHScaleFactor(ScaleFactor):
  '''ScaleFactor class for a spherical harmonic absorption correction'''
  def __init__(self, initial_value, n_parameters, scaling_options=None):
    super(SHScaleFactor, self).__init__(
      initial_value, n_parameters, scaling_options)
    self._harmonic_values = None
    #self.calculate_scales_and_derivatives()

  @property
  def harmonic_values(self):
    return self._harmonic_values

  def update_reflection_data(self, harmonic_values):
    self._harmonic_values = harmonic_values
    self.calculate_scales_and_derivatives()

  #@harmonic_values.setter
  #def harmonic_values(self, values):
  #  '''set a new spherical harmonic coefficient matrix'''
  #  self._harmonic_values = values
  #  self.calculate_scales_and_derivatives()

  def calculate_scales_and_derivatives(self):
    '''calculation of scale factors and derivatives from
       spherical harmonic coefficient matrix'''
    abs_scale = flex.double([1.0] * self._harmonic_values.n_rows)#unity term
    for i, col in enumerate(self._harmonic_values.cols()):
      abs_scale += flex.double(col.as_dense_vector() * self._parameters[i])
    self._inverse_scales = abs_scale
    self._derivatives = self._harmonic_values


class SmoothScaleFactor_2D(SmoothScaleFactor):
  def __init__(self, initial_value, n1_parameters, n2_parameters, scaling_options=None):
    n_parameters = n1_parameters * n2_parameters
    SmoothScaleFactor.__init__(self, initial_value, n_parameters, scaling_options)
    self.n1_parameters = n1_parameters
    self.n2_parameters = n2_parameters
    self.Vr = 0.5
    self.weightsum = None
    self.smoothing_window = 1.5#must be less than 2 to avoid indexing errors

  def set_normalised_values(self, normalised_values1, normalised_values2):
    '''normalised_values is the column from the reflection table
    of the normalised time/resolution etc'''
    self.normalised_values = [normalised_values1, normalised_values2]
    self.scales = flex.double([1.0]*len(self.normalised_values[0]))

  def calculate_smooth_scales(self):
    self.weightsum = flex.double([0.0]*len(self.normalised_values[0]))
    for i, relative_pos_1 in enumerate(self.normalised_values[0]):
      relative_pos_2 = self.normalised_values[1][i]
      max_range_1_to_include = int(relative_pos_1 + self.smoothing_window)
      max_range_2_to_include = int(relative_pos_2 + self.smoothing_window)
      min_range_2_to_include = int((relative_pos_2 - self.smoothing_window)//1) + 1
      min_range_1_to_include = int((relative_pos_1 - self.smoothing_window)//1) + 1
      scale = 0.0
      weightsum = 0.0
      for j in range(min_range_1_to_include, max_range_1_to_include + 1):
        for k in range(min_range_2_to_include, max_range_2_to_include + 1):
          square_distance_to_point = ((float(k) - relative_pos_2)**2
                                      + (float(j) - relative_pos_1)**2)
          if square_distance_to_point < (self.smoothing_window**2):
            scale += (self.scale_factors[(j+1) + ((k+1)*self.n1_parameters)]
                      * np.exp(- square_distance_to_point / self.Vr))
            weightsum += np.exp(-square_distance_to_point/ self.Vr)
      self.weightsum[i] = weightsum
      self.scales[i] = scale/weightsum
    return self.scales

  def calculate_smooth_derivatives(self):
    if not self.weightsum:
      self.calculate_smooth_scales()
    n = len(self.get_normalised_values()[0])
    self.derivatives = flex.double([0.0] * len(self.scale_factors) * n)
    #smoothing_window = 1.5#must be less than 2 to avoid indexing errors
    for i, relative_pos_1 in enumerate(self.normalised_values[0]):
      relative_pos_2 = self.normalised_values[1][i]
      max_range_1_to_include = int(relative_pos_1 + self.smoothing_window)
      max_range_2_to_include = int(relative_pos_2 + self.smoothing_window)
      min_range_2_to_include = int((relative_pos_2 - self.smoothing_window)//1) + 1
      min_range_1_to_include = int((relative_pos_1 - self.smoothing_window)//1) + 1
      for j in range(min_range_1_to_include, max_range_1_to_include + 1):
          for k in range(min_range_2_to_include, max_range_2_to_include + 1):
            square_distance_to_point = ((float(k) - relative_pos_2)**2
                                        + (float(j) - relative_pos_1)**2)
            if square_distance_to_point < (self.smoothing_window**2):
              self.derivatives[(((j+1) + ((k+1)*self.n1_parameters))*n) + i] += (
                np.exp(- square_distance_to_point / self.Vr))/self.weightsum[i]
    return self.derivatives

class SmoothScaleFactor_GridAbsorption(SmoothScaleFactor):
  def __init__(self, initial_value, n1_parameters, n2_parameters, n3_parameters,
               scaling_options=None):
    n_parameters = n1_parameters * n2_parameters * n3_parameters
    SmoothScaleFactor.__init__(self, initial_value, n_parameters, scaling_options)
    self.nx_parameters = n1_parameters
    self.ny_parameters = n2_parameters
    self.ntime_parameters = n3_parameters
    self.Vr = 0.5
    self.weightsum = None

  def set_normalised_values(self, normalised_values1, normalised_values2, normalised_values3):
    '''normalised_values is the column from the reflection table
    of the normalised time/resolution etc'''
    self.normalised_values = [normalised_values1, normalised_values2, normalised_values3]
    self.scales = flex.double([1.0]*len(self.normalised_values[0]))

  def calculate_smooth_scales(self):
    smoothing_window = 1.0#must be 1.0 or less to avoid indexing errors,
    #should probably be fixed or tightly constrained
    self.weightsum = flex.double([0.0]*len(self.normalised_values[0]))
    for datapoint_idx, relative_pos_1 in enumerate(self.normalised_values[0]):
      relative_pos_2 = self.normalised_values[1][datapoint_idx]
      relative_pos_3 = self.normalised_values[2][datapoint_idx]
      max_range_1_to_include = int(relative_pos_1 + smoothing_window)
      max_range_2_to_include = int(relative_pos_2 + smoothing_window)
      max_range_3_to_include = int(relative_pos_3 + smoothing_window)
      min_range_1_to_include = int((relative_pos_1 - smoothing_window)//1) + 1
      min_range_2_to_include = int((relative_pos_2 - smoothing_window)//1) + 1
      min_range_3_to_include = int((relative_pos_3 - smoothing_window)//1) + 1
      scale = 0.0
      weightsum = 0.0
      for i in range(min_range_3_to_include, max_range_3_to_include + 1):
        for j in range(min_range_1_to_include, max_range_1_to_include + 1):
          for k in range(min_range_2_to_include, max_range_2_to_include + 1):
            square_distance_to_point = ((float(k) - relative_pos_2)**2
              + (float(j) - relative_pos_1)**2 + (float(i) - relative_pos_3)**2)
            if square_distance_to_point < (smoothing_window**2):
              scale_idx = (j + (k*self.nx_parameters)) + (i*self.nx_parameters*self.ny_parameters)
              scale += (self.scale_factors[scale_idx] *
                np.exp(- square_distance_to_point / self.Vr))
              weightsum += np.exp(-square_distance_to_point/ self.Vr)
      self.weightsum[datapoint_idx] = weightsum
      self.scales[datapoint_idx] = scale/weightsum
    return self.scales

  def calculate_smooth_derivatives(self):
    if not self.weightsum:
      self.calculate_smooth_scales()
    n = len(self.get_normalised_values()[0])
    self.derivatives = flex.double([0.0] * len(self.scale_factors) * n)
    smoothing_window = 1.0#must be 1.0 or less to avoid indexing problems later
    for datapoint_idx, relative_pos_1 in enumerate(self.normalised_values[0]):
      relative_pos_2 = self.normalised_values[1][datapoint_idx]
      relative_pos_3 = self.normalised_values[2][datapoint_idx]
      max_range_1_to_include = int(relative_pos_1 + smoothing_window)
      max_range_2_to_include = int(relative_pos_2 + smoothing_window)
      max_range_3_to_include = int(relative_pos_3 + smoothing_window)
      min_range_1_to_include = int((relative_pos_1 - smoothing_window)//1) + 1
      min_range_2_to_include = int((relative_pos_2 - smoothing_window)//1) + 1
      min_range_3_to_include = int((relative_pos_3 - smoothing_window)//1) + 1
      for i in range(min_range_3_to_include, max_range_3_to_include + 1):
        for j in range(min_range_1_to_include, max_range_1_to_include + 1):
          for k in range(min_range_2_to_include, max_range_2_to_include + 1):
            square_distance_to_point = ((float(k) - relative_pos_2)**2
              + (float(j) - relative_pos_1)**2 + (float(i) - relative_pos_3)**2)
            if square_distance_to_point < (smoothing_window**2):
              deriv_idx = (j + (k*self.nx_parameters)) + (i*self.nx_parameters*self.ny_parameters)
              self.derivatives[(deriv_idx * n) + datapoint_idx] += (
                np.exp(- square_distance_to_point / self.Vr))/self.weightsum[datapoint_idx]
    return self.derivatives
