#!/usr/bin/env python
#
# dials.compute_delta_cchalf
#
#  Copyright (C) 2018 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#
# LIBTBX_SET_DISPATCHER_NAME dev.dials.compute_delta_cchalf

from __future__ import absolute_import, division, print_function
import sys
from iotbx.reflection_file_reader import any_reflection_file
from iotbx import merging_statistics
from matplotlib import pylab
from matplotlib import cm
from math import sqrt, floor
from cctbx import miller
from cctbx import crystal
from collections import defaultdict


class ResolutionBinner(object):
  '''
  A class to bin the data by resolution

  '''

  def __init__(self, unit_cell, dmin, dmax, nbins):
    '''
    Initialise the binner

    :param unit_cell: A unit cell object
    :param dmin: The maximum resolution
    :param dmax: The minimum resolution
    :param nbins: The number of bins

    '''
    print("Resolution bins")
    assert dmin < dmax
    dmin_inv_sq = 1.0 / dmin**2
    dmax_inv_sq = 1.0 / dmax**2
    self._unit_cell = unit_cell
    self._nbins = nbins
    self._xmin = dmax_inv_sq
    self._xmax = dmin_inv_sq
    assert self._xmin < self._xmax
    self._bin_range = self._xmax - self._xmin
    self._bin_size = self._bin_range / self._nbins
    self._bins = []
    for i in range(self._nbins):
      b0, b1 = self._xmin + i * self._bin_size, self._xmin + (i+1)*self._bin_size
      print("%d: %.3f, %.3f" % (i, sqrt(1/b0**2), sqrt(1/b1**2)))
      self._bins.append((b0,b1))

  def nbins(self):
    '''
    :returns: The number of bins

    '''
    return self._nbins

  def index(self, h):
    '''
    Get the bin index from the miller index

    :param h: The miller index
    :returns: The bin index

    '''
    d = self._unit_cell.d(h)
    d2 = 1/d**2
    bin_index = int(floor((d2 - self._xmin) / self._bin_size))
    if bin_index >= self._nbins:
      bin_index = self._nbins-1
    if bin_index < 0:
      bin_index = 0
    return bin_index


class ReflectionSum(object):
  '''
  A helper class to store sums of X and X**2

  '''
  def __init__(self, sum_x=0, sum_x2=0, n=0):
    self.sum_x = sum_x
    self.sum_x2 = sum_x2
    self.n = n


class BinData(object):
  '''
  A helper class to store mean and variance

  '''
  def __init__(self):
    self.mean = []
    self.var = []


def compute_cchalf(mean, var):
  '''
  Compute the CC 1/2 using the formular from Assmann, Brehm and Diederichs 2016

  :param mean: The list of mean intensities
  :param var: The list of variances on the half set of mean intensities
  :returns: The CC 1/2
  
  '''
  assert len(mean) == len(var)
  n = len(mean)
  mean_of_means = sum(mean) / n
  sigma_e = sum(var) / n
  sigma_y = sum([(m - mean_of_means)**2 for m in mean]) / (n-1)
  cchalf = (sigma_y - sigma_e) / (sigma_y + sigma_e)
  return cchalf


def compute_mean_cchalf_in_bins(bin_data):
  '''
  Compute the mean cchalf averaged across resolution bins
  
  :param bin_data: The mean and variance in each bin
  :returns: The mean CC 1/2

  '''
  mean_cchalf = 0
  count = 0
  for i in range(len(bin_data)):
    mean = bin_data[i].mean
    var = bin_data[i].var
    n = len(mean)
    if n > 1:
      cchalf = compute_cchalf(mean, var)
      mean_cchalf += n*cchalf
      count += n
  mean_cchalf /= count
  return mean_cchalf


class PerImageCChalfStatistics(object):
  '''
  A class to compute per image CC 1/2 statistics

  '''

  def __init__(self, 
               miller_index,
               batch,
               intensity,
               variance,
               unit_cell,
               space_group,
               nbins=10, 
               dmin=None, 
               dmax=None):
    '''
    Initialise

    :param miller_index: The list of miller indices
    :param batch: The list of batch numbers
    :param intensity: The list of intensities
    :param variance: The list of variances
    :param unit_cell: The unit cell or list of unit cells
    :param space_group: The space group
    :param nbins: The number of bins
    :param dmin: The maximum resolution
    :param dmax: The minimum resolution
      
    '''
    
    # Reject reflections with negative variance
    selection = variance > 0
    miller_index = miller_index.select(selection)
    batch = batch.select(selection)
    intensity = intensity.select(selection)
    variance = variance.select(selection)

    # Map the miller indices to the ASU
    cs = crystal.symmetry(unit_cell, space_group=space_group)
    ms = miller.set(cs, miller_index)
    ms_asu = ms.map_to_asu()
    miller_index = ms_asu.indices()

    # Save the arrays
    self._miller_index = miller_index
    self._batch = batch
    self._intensity = intensity
    self._variance = variance

    # Create the resolution bins
    D = [unit_cell.d(h) for h in miller_index]
    self._num_bins = nbins
    self._dmin = dmin
    self._dmax = dmax
    if dmin is None:
      self._dmin = min(D)
    if dmax is None:
      self._dmax = max(D)
    binner = ResolutionBinner(unit_cell, self._dmin, self._dmax, nbins)

    # Create lookups for elements by miller index
    index_lookup = defaultdict(list)
    for i, h in enumerate(miller_index):
      index_lookup[h].append(i)
    
    # Compute the Overall Sum(X) and Sum(X^2) for each unique reflection
    reflection_sums = defaultdict(ReflectionSum)
    for h in index_lookup.keys():
      I = [intensity[i] for i in index_lookup[h]]
      n = len(I)
      sum_x = sum(I)
      sum_x2 = sum(II**2 for II in I)
      reflection_sums[h] = ReflectionSum(sum_x, sum_x2, n)

    # Compute some numbers
    self._num_batches = len(set(batch))
    self._num_reflections = len(miller_index)
    self._num_unique = len(reflection_sums.keys())
    
    print("")
    print("# Batches: ", self._num_batches)
    print("# Reflections: ", self._num_reflections)
    print("# Unique: ", self._num_unique)
   
    # Compute the CC 1/2 for all the data
    self._cchalf_mean = self._compute_cchalf(reflection_sums, binner)
    print("CC 1/2 mean: %.3f" % (100*self._cchalf_mean))
   
    # Compute the CC 1/2 excluding each batch in turn
    self._cchalf = self._compute_cchalf_excluding_each_batch(
      reflection_sums, binner, miller_index, batch, intensity) 

  def _compute_cchalf(self, reflection_sums, binner):
    '''
    Compute the CC 1/2 by computing the CC 1/2 in resolution bins and then
    computing the weighted mean of the binned CC 1/2 values

    '''
    # Compute Mean and variance of reflection intensities
    bin_data = [BinData() for i in range(binner.nbins())]
    for h in reflection_sums.keys():
      sum_x = reflection_sums[h].sum_x
      sum_x2 = reflection_sums[h].sum_x2
      n = reflection_sums[h].n

      if n > 1:
        mean = sum_x / n
        var = (sum_x2 - (sum_x)**2 / n) / (n-1)
        var = var / (0.5*n)
        index = binner.index(h)
        bin_data[index].mean.append(mean)
        bin_data[index].var.append(var)

    # Compute the mean cchalf in resolution bins
    return compute_mean_cchalf_in_bins(bin_data)

  def _compute_cchalf_excluding_each_batch(self, 
                                           reflection_sums, 
                                           binner,
                                           miller_index,
                                           batch,
                                           intensity):
    '''
    Compute the CC 1/2 with an image excluded.

    For each image, update the sums by removing the contribution from the image
    and then compute the CC 1/2 of the remaining data

    '''
   
    # Create a lookup table for each reflection by batch
    batch_lookup = defaultdict(list)
    for i, b in enumerate(batch):
      batch_lookup[b].append(i)

    # Compute CC1/2 minus each batch
    cchalf_i = {}
    for batch in batch_lookup.keys():

      # Find all observations from this batch and create a lookup based on
      # miller index
      index_lookup = defaultdict(list)
      for i in batch_lookup[batch]:
        index_lookup[miller_index[i]].append(i)

      # Loop through all the reflections and remove the contribution of
      # reflections from the current batch
      batch_reflection_sums = defaultdict(ReflectionSum)
      for h in reflection_sums.keys():
        sum_x = reflection_sums[h].sum_x
        sum_x2 = reflection_sums[h].sum_x2
        n = reflection_sums[h].n
        for i in index_lookup[h]:
          sum_x -= self._intensity[i]
          sum_x2 -= self._intensity[i]**2
          n -= 1
        batch_reflection_sums[h] = ReflectionSum(sum_x, sum_x2, n)
    
      # Compute the CC 1/2 without the reflections from the current batch
      cchalf = self._compute_cchalf(batch_reflection_sums, binner)
      cchalf_i[batch] = cchalf
      print("CC 1/2 excluding batch %d: %.3f" % (batch, 100*cchalf))

    return cchalf_i

  def num_batches(self):
    '''
    Return the number of batches

    '''
    return len(self._cchalf)

  def num_reflections(self):
    '''
    Return the number of reflections

    '''
    return self._num_reflections

  def num_unique(self):
    '''
    Return the number of unique reflections

    '''
    return self._num_unique

  def mean_cchalf(self):
    '''
    Return the mean CC 1/2

    '''
    return self._cchalf_mean

  def cchalf_i(self):
    '''
    Return the CC 1/2 for each image excluded

    '''
    return self._cchalf

  def delta_cchalf_i(self):
    '''
    Return the Delta CC 1/2 for each image excluded

    '''
    return dict((k, self._cchalf_mean - v) for k, v in self._cchalf.iteritems())


if __name__ == '__main__':

  # Read the mtz file
  reader = any_reflection_file(sys.argv[1])

  # Get the columns as miller arrays
  miller_arrays = reader.as_miller_arrays(merge_equivalents=False)

  # Select the desired columns
  intensities = None
  batches = None
  for array in miller_arrays:
    if array.info().labels == ['I', 'SIGI']:
      intensities = array
    if array.info().labels == ['BATCH']:
      batches = array
  assert intensities is not None
  assert batches is not None
  assert len(batches.data()) == len(intensities.data())

  # Compute the resolution of each reflection
  unit_cell = intensities.unit_cell()
    
  # Create the statistics object
  statistics = PerImageCChalfStatistics(
    intensities.indices(),
    batches.data(),
    intensities.data(),
    intensities.sigmas()**2,
    unit_cell,
    intensities.crystal_symmetry().space_group())

  # Print out the Batches in order of delta cc 1/2
  delta_cchalf_i = statistics.delta_cchalf_i()
  batches = list(delta_cchalf_i.keys())
  sorted_index = sorted(range(len(batches)), key=lambda x: delta_cchalf_i[batches[x]])
  for i in sorted_index:
    print("Batch: %d, Delta CC 1/2: %.3f" % (batches[i], 100*delta_cchalf_i[batches[i]]))

  # Make a plot of delta cc 1/2
  from matplotlib import pylab
  pylab.hist(delta_cchalf_i.values())
  pylab.xlabel("Delta CC 1/2")
  pylab.show()

  X = list(delta_cchalf_i.keys())
  Y = list(delta_cchalf_i.values())
  pylab.plot(X, Y)
  pylab.xlabel("Batch number")
  pylab.ylabel("Delta CC 1/2")
  pylab.show()
    
