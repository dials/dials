from dials.array_family import flex
import numpy as np
import math as math

'''
These classes define objects to hold a flex array of scale factors,
including the required methods for calculating smooth scale factors and
'''

class ScaleFactor(object):
  def __init__(self, initial_value, n_parameters, scaling_options=None):
    self.scale_factors = flex.double([initial_value] * n_parameters)
    self.n_params = n_parameters

  def set_scale_factors(self, scale_factors):
    if len(scale_factors) != len(self.scale_factors):
      assert 0, '''attempting to set a new set of scale factors of different
      length than previous assignment: was %s, attempting %s''' % (
        len(self.scale_factors), len(scale_factors))
    self.scale_factors = scale_factors

  def get_scale_factors(self):
    return self.scale_factors

class K_ScaleFactor(ScaleFactor):
  def __init__(self, initial_value, n_refl, scaling_options=None):
    ScaleFactor.__init__(self, initial_value, 1, scaling_options)
    self.n_refl = n_refl

  def set_n_refl(self, n_refl):
    self.n_refl = n_refl

  def get_scales_of_reflections(self):
    '''only one parameter for a K_ScaleFactor'''
    return flex.double([self.scale_factors[0]]*self.n_refl)

  def get_derivatives(self):
    'return derivatives'
    return flex.double([1.0]*self.n_refl)

class B_ScaleFactor(ScaleFactor):
  def __init__(self, initial_value, d_values, scaling_options=None):
    ScaleFactor.__init__(self, initial_value, 1, scaling_options)
    self.d_values = d_values

  def set_d_values(self, d_values):
    self.d_values = d_values

  def get_scales_of_reflections(self):
    '''only one parameter for a B_ScaleFactor'''
    scales = flex.double([self.scale_factors[0]]*len(self.d_values))
    return flex.double(np.exp(scales/(2.0 * (self.d_values**2))))

  def get_derivatives(self):
    'return derivatives'
    scales = self.get_scales_of_reflections()
    return scales / (2.0 * (self.d_values**2))
 

class SmoothScaleFactor(ScaleFactor):
  def __init__(self, initial_value, n_parameters, scaling_options=None):
    ScaleFactor.__init__(self, initial_value, n_parameters, scaling_options)
    self.normalised_values = None
    self.Vr = 1.0# 0.7

  def set_normalised_values(self, normalised_values):
    '''normalised_values is the column from the reflection table
    of the normalised time/resolution etc'''
    self.normalised_values = normalised_values
    self.scales = flex.double([1.0]*len(self.normalised_values))

  def get_normalised_values(self):
    return self.normalised_values

  def get_scales_of_reflections(self):
    return self.scales

class SmoothScaleFactor_1D(SmoothScaleFactor):
  def __init__(self, initial_value, n_parameters, scaling_options=None):
    SmoothScaleFactor.__init__(self, initial_value, n_parameters, scaling_options)
    self.smoothing_window = 2.5 #must be less than 3 to avoid indexing errors

  def set_normalised_values(self, normalised_values):
    '''normalised_values is the column from the reflection table
    of the normalised time/resolution etc'''
    self.normalised_values = normalised_values
    self.n_norm_val = int(len(self.normalised_values))
    self.max_scale_to_include = flex.int([])
    self.min_scale_to_include = flex.int([])
    self.scales = flex.double([1.0]*self.n_norm_val)
    for norm_val in self.normalised_values:
      self.max_scale_to_include.append(int(norm_val + self.smoothing_window))
      self.min_scale_to_include.append(int((norm_val - self.smoothing_window)//1) + 1)
    self.zip_vals = zip(self.normalised_values, 
      self.min_scale_to_include, self.max_scale_to_include)
    self.exponentials = flex.float([])
    for norm_val, min_s, max_s in self.zip_vals:
      for j in range(min_s, max_s + 1):
        self.exponentials.append(math.exp(-((norm_val - float(j))**2) / self.Vr))

  def calculate_smooth_scales(self):
    self.weightsum = flex.float([])
    self.scales = flex.float([])
    counter = 0
    for norm_val, min_s, max_s in self.zip_vals:
      scale = 0.0
      weightsum = 0.0
      for j in range(min_s, max_s + 1):
        scale += self.scale_factors[j+2] * self.exponentials[counter]
        weightsum += self.exponentials[counter]
        counter += 1
      self.weightsum.append(weightsum)
      self.scales.append(scale/weightsum)
    return self.scales.as_double()

  def calculate_smooth_derivatives(self):
    self.derivatives = flex.float([0.0] * len(self.scale_factors) * self.n_norm_val)
    counter = 0
    for i, (norm_val, min_s, max_s) in enumerate(self.zip_vals):
      for j in range(min_s, max_s + 1):
        self.derivatives[((j+2)*self.n_norm_val)+i] += (
          self.exponentials[counter]/self.weightsum[i])
        counter += 1
    return self.derivatives.as_double()

class SmoothScaleFactor_1D_Bfactor(SmoothScaleFactor_1D):
  def __init__(self, initial_value, n_parameters, d_values, scaling_options=None):
    SmoothScaleFactor_1D.__init__(self, initial_value, n_parameters, scaling_options)
    self.d_values = d_values
    self.Vr = 0.5

  def set_d_values(self, d_values):
    if len(d_values) != len(self.normalised_values):
      assert 0, '''attempting to set a new set of d factors of different
      length to the data: %s vs %s''' % (
        len(d_values), len(self.normalised_values))
    self.d_values = d_values

  def calculate_smooth_scales(self):
    self.weightsum = flex.float([])
    self.scales = flex.float([])
    counter = 0
    for norm_val, min_s, max_s in self.zip_vals:
      scale = 0.0
      weightsum = 0.0
      for j in range(min_s, max_s + 1):
        scale += self.scale_factors[j+2] * self.exponentials[counter]
        weightsum += self.exponentials[counter]
        counter += 1
      self.weightsum.append(weightsum)
      self.scales.append(scale/weightsum)
    self.scales = flex.double(np.exp(self.scales.as_double()/(2.0 * (self.d_values**2))))
    return self.scales

  def calculate_smooth_derivatives(self):
    self.derivatives = flex.float([0.0] * len(self.scale_factors) * self.n_norm_val)
    counter = 0
    for i, (norm_val, min_s, max_s) in enumerate(self.zip_vals):
      for j in range(min_s, max_s + 1):
        self.derivatives[((j+2)*self.n_norm_val)+i] += (
          self.exponentials[counter]/self.weightsum[i])
        counter += 1
    multiplicative_factor = (flex.double(np.exp(self.scales.as_double()/(2.0 * (self.d_values**2))))
                             / (2.0 * (self.d_values**2))) 
    return self.derivatives.as_double() * flex.double(np.tile(multiplicative_factor, self.n_params))


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
    smoothing_window = 1.5#must be less than 2 to avoid indexing errors
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

class SphericalAbsorption_ScaleFactor(ScaleFactor):
  def __init__(self, initial_value, n_parameters, values, scaling_options=None):
    ScaleFactor.__init__(self, initial_value, n_parameters, scaling_options)
    self.harmonic_values = values
    self.scales = self.calculate_smooth_scales()
    self.derivatives = self.calculate_smooth_derivatives()

  def set_values(self, values):
    self.harmonic_values = values
    self.scales = self.calculate_smooth_scales()
    self.derivatives = self.calculate_smooth_derivatives()

  def get_values(self):
    return self.harmonic_values

  def get_scales_of_reflections(self):
    return self.scales

  def calculate_smooth_scales(self):
    '''name the methods the same as the rest of the aimless scale factors
    to allow generic calling - should all the names just be the same, no smooth?'''
    abs_scale = flex.double([1.0]*len(self.harmonic_values))
    for n in range(self.harmonic_values.ncols()):
      abs_scale += self.harmonic_values[str(n)] * self.scale_factors[n]
    return abs_scale

  def calculate_smooth_derivatives(self):
    derivatives = flex.double([])
    for n in range(self.harmonic_values.ncols()):
      derivatives.extend(self.harmonic_values[str(n)])
    return derivatives
