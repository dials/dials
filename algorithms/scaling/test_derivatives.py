from cctbx.array_family import flex
import numpy as np
from scitbx.random import variate, poisson_distribution
from libtbx.test_utils import approx_equal


class test_dataset(object):
  '''simple empty 'data-manager-like' object on which
  gradients can be calculated'''
  def __init__(self):
    self.intensities = None
    self.scale_factors = None
    self.Ih_values = None
    self.Ih_array = None
    self.h_index_counter_array = None
    self.h_index_cumulative_array = None
    self.scaleweights = None
    self.bin_indices = None
    self.n_params = None
    self.derivatives = None

  def calc_Ih(self):
    '''method to calculate the best estimate for the intensity for
    unique reflections'''
    gsq = (((self.scale_factors)**2) * self.scaleweights)
    sumgsq = flex.double(np.add.reduceat(gsq, self.h_index_cumulative_array[:-1]))
    gI = ((self.scale_factors * self.intensities) * self.scaleweights)
    sumgI = flex.double(np.add.reduceat(gI, self.h_index_cumulative_array[:-1]))
    sumweights = flex.double(np.add.reduceat(self.scaleweights, self.h_index_cumulative_array[:-1]))
    self.Ih_array = flex.double([val/ sumgsq[i] if sumweights[i] > 0.0
                                 else 0.0 for i, val in enumerate(sumgI)])
    self.Ih_values = flex.double(np.repeat(self.Ih_array, self.h_index_counter_array))

  def calculate_gradient_fd(self):
    '''calculate gradient array with finite difference approach'''
    delta = 1.0e-6
    gradients = [0.0] * self.n_params
    #iterate over number of parameters and calculate the gradient for each
    for i in range(self.n_params):
      intensities = self.intensities
      weights = self.scaleweights

      #calculate R_low
      for j in range(len(intensities)):
        if self.bin_indices[j] == i:
          self.scale_factors[j] -= (0.5 * delta)
      self.calc_Ih()
      Ih_values_low = self.Ih_values
      scale_factors_low = self.scale_factors
      R_low = ((((intensities - (scale_factors_low * Ih_values_low))**2) * weights))

      #calculate R_upper
      for j in range(len(intensities)):
        if self.bin_indices[j] == i:
          self.scale_factors[j] += (delta)
      self.calc_Ih()
      Ih_values_upper = self.Ih_values
      scale_factors_upper = self.scale_factors
      R_upper = ((((intensities - (scale_factors_upper * Ih_values_upper))**2) * weights))

      #add delta/2 back to scale factors to return them to their initial values
      for j in range(len(intensities)):
        if self.bin_indices[j] == i:
          self.scale_factors[j] -= (0.5 * delta)
      self.calc_Ih()

      gradients[i] += (flex.sum(R_upper) - flex.sum(R_low)) / delta
    return gradients

  def calculate_gradient(self):
    '''calculate gradient array using the expression derived from differentiation
    of the least-squares target'''
    gradient = []
    gsq = ((self.scale_factors)**2) * self.scaleweights
    sumgsq = np.add.reduceat(gsq, self.h_index_cumulative_array[:-1])
    rhl = self.intensities - (self.Ih_values * self.scale_factors)
    num = len(self.intensities)
    dIh = ((self.intensities * self.scaleweights)
           - (self.Ih_values * 2.0 * self.scale_factors * self.scaleweights))

    for i in range(self.n_params):
      dIh_g = (dIh * self.derivatives[i*num:(i+1)*num])
      dIh_g = np.add.reduceat(dIh_g, self.h_index_cumulative_array[:-1])/sumgsq
      dIh_g = flex.double(np.repeat(dIh_g, self.h_index_counter_array))
      drdp = -((self.Ih_values * self.derivatives[i*num:(i+1)*num])
               + (self.scale_factors * dIh_g))
      grad = (2.0 * rhl * self.scaleweights * drdp)
      gradient.append(flex.sum(grad))
    return gradient


if __name__ == '__main__':
  data = test_dataset()
  data.n_params = 4
  n_obs = 10

  #generate data
  K = [1 + e/40. for e in range(n_obs)]
  data.scale_factors = flex.double([1./e for e in K])
  means = [100 * e for e in K]
  data.intensities = flex.double([variate(poisson_distribution(e))() for e in means])
  data.scaleweights = flex.double([1./e for e in means])

  #simulate a parameterisation - bin randomly into four different bins
  data.bin_indices = flex.int([0, 3, 1, 3, 2, 0, 2, 0, 2, 3])
  #simulate the grouping into groups of unique reflections
  data.h_index_cumulative_array = flex.double([0, 2, 5, 7, 10])
  data.h_index_counter_array = flex.double([2, 3, 2, 3])

  #initialise derivates of the scales w.r.t parameters for the full gradient
  #calculation. In its current form this is the xds-like parameterisation.
  data.derivatives = flex.double([0.0] * n_obs * data.n_params)
  print list(data.derivatives)
  for j, a in enumerate(data.bin_indices):
    data.derivatives[(a*n_obs)+j] = 1.0

  #needs an initial value of Ih for calculating gradient
  data.calc_Ih()

  #calculate the gradients by direct calculation and finite difference method
  print data.calculate_gradient()
  print data.calculate_gradient_fd()
