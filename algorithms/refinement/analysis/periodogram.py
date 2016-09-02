#!/usr/bin/env cctbx.python
#
#  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#
"""Calculate the periodogram of real evenly-spaced data"""

from __future__ import division
from dials.array_family import flex
from scitbx import fftpack

class Periodogram(object):
  """Calculate the raw periodogram of real evenly-spaced data. This class gives
  the same spectrum as R's spec.pgram function when called using
  spec.pgram(x, detrend=F, taper=0, fast=F)"""

  def __init__(self, x, demean=True, detrend=True):

    # Ensure x is copied as it will be changed in-place
    x = flex.double(x).deep_copy()
    n = len(x)

    if detrend:
      t = flex.size_t_range(n).as_double() + 1 - (n + 1)/2
      inv_sumt2 = 1./t.dot(t)
      x = x - flex.mean(x) - x.dot(t) * t * inv_sumt2
    elif demean:
      x -= flex.mean(x)

    # determine frequencies
    stop = ((n - (n % 2)) // 2) + 1
    self.freq = flex.double([i / n for i in range(1, stop)])

    fft = fftpack.real_to_complex(n)
    n = fft.n_real()
    m = fft.m_real()
    x.extend(flex.double(m-n, 0.))
    xf = fft.forward(x)

    # remove the first term (DC offset), get abs length of complex number and
    # normalise by n to produce the raw periodogram
    self.spec = flex.norm(xf[1:]) / n

  def plot(self, show=True):

    import matplotlib.pyplot as plt
    line, = plt.semilogy(self.freq, self.spec)
    plt.xlabel('frequency')
    plt.ylabel('spectrum')
    if show: plt.show()
    return plt

