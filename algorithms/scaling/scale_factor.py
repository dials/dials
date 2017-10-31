from dials.array_family import flex
import numpy as np

'''
These classes define objects to hold a flex array of scale factors,
including the required methods for calculating smooth scale factors and
'''

class ScaleFactor(object):
  def __init__(self, initial_value, n_parameters, scaling_options=None):
    self.scale_factors = flex.double([initial_value] * n_parameters)

  def set_scale_factors(self, scale_factors):
    if len(scale_factors) != len(self.scale_factors):
      assert 0, '''attempting to set a new set of scale factors of different
      length than previous assignment: was %s, attempting %s''' % (
        len(self.scale_factors), len(scale_factors))
    self.scale_factors = scale_factors

  def get_scale_factors(self):
    return self.scale_factors

class B_ScaleFactor(ScaleFactor):
  def __init__(self, initial_value, n_parameters, d_values, scaling_options=None):
    ScaleFactor.__init__(self, initial_value, n_parameters, scaling_options)
    self.d_values = d_values

  def set_d_values(self, d_values):
    self.d_values = d_values
  
      

class SmoothScaleFactor(ScaleFactor):
  def __init__(self, initial_value, n_parameters, scaling_options=None):
    ScaleFactor.__init__(self, initial_value, n_parameters, scaling_options)
    self.normalised_values = None
    self.Vr = 0.7
    self.problim = 4.0
    if scaling_options:
      self.Vr = scaling_options['Vr']
      self.problim = scaling_options['problim']

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
    if scaling_options:
      self.Vr = scaling_options['Vr']
      self.problim = scaling_options['problim']

  def calculate_smooth_scales(self):
    (Vr, problim) = (self.Vr, self.problim)
    smoothing_window = 2.5#must be less than 3 to avoid indexing errors
    self.weightsum = flex.double([0.0]*len(self.normalised_values))
    for i, relative_pos in enumerate(self.normalised_values):
      max_scale_to_include = int(relative_pos + smoothing_window)
      min_scale_to_include = int((relative_pos - smoothing_window)//1) + 1
      scale = 0.0
      weightsum = 0.0
      for j in range(min_scale_to_include, max_scale_to_include + 1):
        scale += self.scale_factors[j+2] * np.exp(-((relative_pos - float(j))**2) / Vr)
        weightsum += np.exp(-((relative_pos - float(j))**2) / Vr)
      self.weightsum[i] = weightsum
      self.scales[i] = scale/weightsum
    return self.scales

  def calculate_smooth_derivatives(self):
    n = len(self.normalised_values)
    self.derivatives = flex.double([0.0] * len(self.scale_factors) * n)
    (Vr, problim) = (self.Vr, self.problim)
    smoothing_window = 2.5#must be less than 3 to avoid indexing errors
    for i, relative_pos in enumerate(self.normalised_values):
      max_scale_to_include = int(relative_pos + smoothing_window)
      min_scale_to_include = int((relative_pos - smoothing_window)//1) + 1
      for j in range(min_scale_to_include, max_scale_to_include + 1):
        self.derivatives[((j+2)*n)+i] += (#self.scale_factors[j+1] *
          np.exp(-((relative_pos - float(j))**2) / Vr))/self.weightsum[i]
    return self.derivatives

class SmoothScaleFactor_1D_Bfactor(SmoothScaleFactor_1D):
  def __init__(self, initial_value, n_parameters, d_values, scaling_options=None):
    SmoothScaleFactor_1D.__init__(self, initial_value, n_parameters, scaling_options)
    self.d_values = d_values

  def set_d_values(self, d_values):
    if len(d_values) != len(self.normalised_values):
      assert 0, '''attempting to set a new set of d factors of different
      length to the data: %s vs %s''' % (
        len(d_values), len(self.normalised_values))
    self.d_values = d_values

  def calculate_smooth_scales(self):
    (Vr, problim) = (self.Vr, self.problim)
    smoothing_window = 2.5#must be less than 3 to avoid indexing errors
    self.scales = flex.double([0.0]*len(self.normalised_values))
    self.weightsum = flex.double([0.0]*len(self.normalised_values))
    for i, relative_pos in enumerate(self.normalised_values):
      max_scale_to_include = int(relative_pos + smoothing_window)
      min_scale_to_include = int((relative_pos - smoothing_window)//1) + 1
      scale = 0.0
      weightsum = 0.0
      for j in range(min_scale_to_include, max_scale_to_include + 1):
        scale += self.scale_factors[j+2] * np.exp(-((relative_pos - float(j))**2) / Vr)
        weightsum += np.exp(-((relative_pos - float(j))**2) / Vr)
      self.weightsum[i] = weightsum
      self.scales[i] = scale/weightsum
    self.scales = flex.double(np.exp(self.scales/(2.0 * (self.d_values**2))))
    return self.scales

class SmoothScaleFactor_2D(SmoothScaleFactor):
  def __init__(self, initial_value, n1_parameters, n2_parameters, scaling_options=None):
    n_parameters = n1_parameters * n2_parameters
    SmoothScaleFactor.__init__(self, initial_value, n_parameters, scaling_options)
    self.n1_parameters = n1_parameters
    self.n2_parameters = n2_parameters
    if scaling_options:
      self.Vr = scaling_options['Vr']
      self.problim = scaling_options['problim']
    self.Vr = 0.5
    self.weightsum = None

  def set_normalised_values(self, normalised_values1, normalised_values2):
    '''normalised_values is the column from the reflection table
    of the normalised time/resolution etc'''
    self.normalised_values = [normalised_values1, normalised_values2]
    self.scales = flex.double([1.0]*len(self.normalised_values[0]))

  def calculate_smooth_scales(self):
    (Vr, problim) = (self.Vr, self.problim)
    smoothing_window = 1.5#must be less than 2 to avoid indexing errors
    self.weightsum = flex.double([0.0]*len(self.normalised_values[0]))
    for i, relative_pos_1 in enumerate(self.normalised_values[0]):
      relative_pos_2 = self.normalised_values[1][i]
      max_range_1_to_include = int(relative_pos_1 + smoothing_window)
      max_range_2_to_include = int(relative_pos_2 + smoothing_window)
      min_range_2_to_include = int((relative_pos_2 - smoothing_window)//1) + 1
      min_range_1_to_include = int((relative_pos_1 - smoothing_window)//1) + 1
      scale = 0.0
      weightsum = 0.0
      for j in range(min_range_1_to_include, max_range_1_to_include + 1):
        for k in range(min_range_2_to_include, max_range_2_to_include + 1):
          square_distance_to_point = ((float(k) - relative_pos_2)**2
                                      + (float(j) - relative_pos_1)**2)
          if square_distance_to_point < (smoothing_window**2):
            scale += (self.scale_factors[(j+1) + ((k+1)*self.n1_parameters)]
                      * np.exp(- square_distance_to_point / Vr))
            weightsum += np.exp(-square_distance_to_point/ Vr)
      self.weightsum[i] = weightsum
      self.scales[i] = scale/weightsum
    return self.scales

  def calculate_smooth_derivatives(self):
    if not self.weightsum:
      self.calculate_smooth_scales()
    n = len(self.get_normalised_values()[0])
    self.derivatives = flex.double([0.0] * len(self.scale_factors) * n)
    (Vr, problim) = (self.Vr, self.problim)
    smoothing_window = 1.5#must be less than 2 to avoid indexing errors
    for i, relative_pos_1 in enumerate(self.normalised_values[0]):
      relative_pos_2 = self.normalised_values[1][i]
      max_range_1_to_include = int(relative_pos_1 + smoothing_window)
      max_range_2_to_include = int(relative_pos_2 + smoothing_window)
      min_range_2_to_include = int((relative_pos_2 - smoothing_window)//1) + 1
      min_range_1_to_include = int((relative_pos_1 - smoothing_window)//1) + 1
      for j in range(min_range_1_to_include, max_range_1_to_include + 1):
          for k in range(min_range_2_to_include, max_range_2_to_include + 1):
            square_distance_to_point = ((float(k) - relative_pos_2)**2
                                        + (float(j) - relative_pos_1)**2)
            if square_distance_to_point < (smoothing_window**2):
              self.derivatives[(((j+1) + ((k+1)*self.n1_parameters))*n) + i] += (
                np.exp(- square_distance_to_point / Vr))/self.weightsum[i]
    return self.derivatives

class SmoothScaleFactor_GridAbsorption(SmoothScaleFactor):
  def __init__(self, initial_value, n1_parameters, n2_parameters, n3_parameters, 
               scaling_options=None):
    n_parameters = n1_parameters * n2_parameters * n3_parameters
    SmoothScaleFactor.__init__(self, initial_value, n_parameters, scaling_options)
    self.nx_parameters = n1_parameters
    self.ny_parameters = n2_parameters
    self.ntime_parameters = n3_parameters
    if scaling_options:
      self.Vr = scaling_options['Vr']
      self.problim = scaling_options['problim']
    self.Vr = 0.5
    self.weightsum = None

  def set_normalised_values(self, normalised_values1, normalised_values2, normalised_values3):
    '''normalised_values is the column from the reflection table
    of the normalised time/resolution etc'''
    self.normalised_values = [normalised_values1, normalised_values2, normalised_values3]
    self.scales = flex.double([1.0]*len(self.normalised_values[0]))

  def calculate_smooth_scales(self):
    (Vr, problim) = (self.Vr, self.problim)
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
                np.exp(- square_distance_to_point / Vr))
              weightsum += np.exp(-square_distance_to_point/ Vr)
      self.weightsum[datapoint_idx] = weightsum
      self.scales[datapoint_idx] = scale/weightsum
    return self.scales


  def calculate_smooth_derivatives(self):
    if not self.weightsum:
      self.calculate_smooth_scales()
    n = len(self.get_normalised_values()[0])
    self.derivatives = flex.double([0.0] * len(self.scale_factors) * n)
    (Vr, problim) = (self.Vr, self.problim)
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
                np.exp(- square_distance_to_point / Vr))/self.weightsum[datapoint_idx]
    return self.derivatives

class SphericalAbsorption_ScaleFactor(ScaleFactor):
  def __init__(self, initial_value, n_parameters, values, scaling_options=None):
    ScaleFactor.__init__(self, initial_value, n_parameters, scaling_options)
    self.harmonic_values = values
    self.scales = self.calculate_scales()
    self.derivatives = self.calculate_derivatives()

  def set_values(self, values):
    self.harmonic_values = values
    self.scales = self.calculate_scales()
    self.derivatives = self.calculate_derivatives()

  def get_values(self):
    return self.harmonic_values

  def get_scales_of_reflections(self):
    return self.scales

  def calculate_scales(self):
    abs_scale = flex.double([1.0]*len(self.harmonic_values))
    count = 0
    for key in self.harmonic_values.keys():
      abs_scale += self.harmonic_values[key] * self.scale_factors[count]
      count += 1
    return abs_scale

  def calculate_derivatives(self):
    derivatives = flex.double([])
    for key in self.harmonic_values.keys():
      derivatives.extend(self.harmonic_values[key])
    return derivatives
