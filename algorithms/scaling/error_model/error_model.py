"""
Error model classes for scaling.
"""
import logging
from dials.array_family import flex
from scitbx import sparse
logger = logging.getLogger('dials')

class BasicErrorModel(object):
  """
  Object to manage calculation of deviations for an error model.
  """

  def __init__(self, Ih_table):
    logger.info("Initialising an error model for refinement.")
    self.Ih_table = Ih_table
    self.n_h = self.Ih_table.calc_nh()
    self.sigmaprime = None
    self.delta_hl = None
    self.bin_variances = None
    self._summation_matrix = self.create_summation_matrix()
    self._bin_counts = flex.double(self.Ih_table.size, 1.0) * self.summation_matrix
    self.refined_parameters = None

  @property
  def summation_matrix(self):
    """A sparse matrix to allow summation over intensity groups."""
    return self._summation_matrix

  @property
  def bin_counts(self):
    """An array of the number of intensities assigned to each bin."""
    return self._bin_counts

  def calc_sigmaprime(self, x):
    """Calculate the error from the model."""
    sigmaprime = x[0] * ((self.Ih_table.variances)
      + ((x[1]*self.Ih_table.intensities)**2))**0.5
    return sigmaprime

  def calc_deltahl(self):
    """Calculate the normalised deviations from the model."""
    I_hl = self.Ih_table.intensities
    g_hl = self.Ih_table.inverse_scale_factors
    I_h = self.Ih_table.Ih_values
    prefactor = ((self.n_h - flex.double(self.n_h.size(), 1.0)) / self.n_h)**0.5
    delta_hl = prefactor * ((I_hl/g_hl) - I_h) / self.sigmaprime
    return delta_hl

  def update_for_minimisation(self, x):
    """"Calculate the updated quantites."""
    self.sigmaprime = self.calc_sigmaprime(x)
    self.delta_hl = self.calc_deltahl()
    self.bin_variances = self.calculate_bin_variances()

  def create_summation_matrix(self):
    """"Create a summation matrix to allow sums into intensity bins."""
    sel = flex.sort_permutation(self.Ih_table.intensities)
    #sel is the list of indices in order of small to high
    n = self.Ih_table.size
    if n < 10000: # what is a sensible limit here?
      n_bins = 10 # what is a sensible number of bins?
    else:
      n_bins = 20
    summation_matrix = sparse.matrix(n, n_bins)
    n_cumul_array = flex.int([])
    for i in range(0, n_bins+1):
      n_cumul_array.append((i*n)//n_bins)
    for j in range(len(n_cumul_array)-1):
      for index in sel[n_cumul_array[j]:n_cumul_array[j+1]]:
        summation_matrix[index, j] = 1
    return summation_matrix

  def calculate_bin_variances(self):
    """Calculate the variance of each bin."""
    sum_deltasq = (self.delta_hl**2) * self.summation_matrix
    sum_delta_sq = (self.delta_hl * self.summation_matrix)**2
    return (sum_deltasq/self.bin_counts) - (sum_delta_sq/(self.bin_counts**2))
