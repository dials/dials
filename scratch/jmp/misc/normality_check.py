from __future__ import division


def Simple(data, n_sigma=3):
  import numpy

  # Get mean/sdev
  mean = numpy.mean(data)
  sdev = numpy.std(data)

  # Get min/max
  mind = numpy.min(data)
  maxd = numpy.max(data)

  # Calculate t-statistic of min/max
  min_n_sigma = (mean - mind) / sdev
  max_n_sigma = (maxd - mean) / sdev

  # return whether within required sigma
  return max_n_sigma < n_sigma and min_n_sigma < n_sigma


if __name__ == '__main__':
  data = [10, 10, 10, 10, 10, 10, 10, 10, 10, 100]
  print Simple(data)
