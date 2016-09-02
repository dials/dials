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

class Kernel(object):
  """A discrete symmetric normalized smoothing kernel for use with kernapply.
  Currently allowed types are 'daniell' or 'modified.daniell'"""

  def __init__(self, name='daniell', m=2):
    if not name in ['daniell', 'modified.daniell']:
      raise TypeError('Unknown kernel type "{0}"'.format(name))

    self.m = m
    self.name = name
    if name == 'daniell': self._set_daniell()
    if name == 'modified.daniell': self._set_modified_daniell()
    return

  def _set_daniell(self):
    self.coef = flex.double(self.m + 1, 1./(2. * self.m + 1))

  def _set_modified_daniell(self):
    self.coef = flex.double(self.m + 1, 1./(2. * self.m))
    self.coef[self.m] = self.coef[self.m] * 0.5

def kernapply(x, k, circular=False):
  """Convolve a sequence x with a Kernel k"""

  x = flex.double(x).deep_copy()
  lenx = len(x)
  w = flex.double(lenx, 0.0)
  w.set_selected(flex.size_t_range(k.m + 1), k.coef)
  sel = lenx -1 - flex.size_t_range(k.m)
  w.set_selected(sel, k.coef[1:])

  # do convolution in the Fourier domain
  fft = fftpack.real_to_complex(lenx)
  n = fft.n_real()
  m = fft.m_real()
  x.extend(flex.double(m-n, 0.))
  w.extend(flex.double(m-n, 0.))
  conv = fft.forward(x) * fft.forward(w)

  # extend result by the reverse conjugate, omitting the DC offset and Nyquist
  # frequencies. If fft.n_real() is odd there is no Nyquist term.
  end = fft.n_complex() - (fft.n_real() + 1) % 2
  conv.extend(flex.conj(conv[1:end]).reversed())

  # transform back, take real part and scale
  fft = fftpack.complex_to_complex(len(conv))
  result = fft.backward(conv).parts()[0] / n

  if circular:
    return result
  else:
    return result[(k.m):(lenx-k.m)]


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

