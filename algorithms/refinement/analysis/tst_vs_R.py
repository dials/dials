#!/usr/bin/env cctbx.python
#
#  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#
"""Compare Periodogram with R's spec.pgram"""

from __future__ import division
try:
  import rpy2.robjects as robjects
except ImportError as e:
  print "This script requires the rpy2 package to be installed."
  raise
from dials.array_family import flex
from scitbx.math.periodogram import Periodogram
from libtbx.test_utils import approx_equal

class rpgram(object):
  """Calculate a raw periodogram using R"""
  def __init__(self, x, detrend=True, spans=None):

    if spans is None: spans = robjects.r("NULL")
    dat = robjects.FloatVector(list(x))
    pgram = robjects.r['spec.pgram']
    result = pgram(dat, spans=spans, detrend=detrend, taper=0, fast=False, plot=False)
    self.result = result
    self.spec = flex.double(result.rx2('spec'))
    self.freq = flex.double(result.rx2('freq'))
    return

  def plot(self):
    rplot = robjects.r('plot')
    rplot(self.result, main="Periodogram")


def test1():
  # Compare over a range of lengths from 2 to 200 with random data
  for i in range(2, 201):
    dat = flex.random_double(i)
    a = Periodogram(dat)
    b = rpgram(dat)

    assert approx_equal(a.freq, b.freq)
    assert approx_equal(a.spec, b.spec)
  print "OK"

def test2():
  # compare plots
  dat = flex.random_double(50)
  a = Periodogram(dat)
  b = rpgram(dat)

  from matplotlib.pyplot import ion
  ion()
  a.plot()
  b.plot()

def test3():
  # compare smoothed pgrams with even and odd length sequences
  for i in range(2):
    dat = flex.random_double(50+i)
    a = Periodogram(dat, spans=4)
    b = rpgram(dat, spans=4)

    #from matplotlib.pyplot import ion
    #ion()
    #a.plot()
    #b.plot()

    assert approx_equal(a.freq, b.freq)
    assert approx_equal(a.spec, b.spec)

    print "OK"

def test4():
  # compare kernapply
  from scitbx.math.periodogram import Kernel, kernapply
  spans = 4
  dat = flex.random_double(50)
  k1 = Kernel('modified.daniell', spans//2)
  a = kernapply(dat, k1, circular=True)

  x = robjects.FloatVector(list(dat))
  kernel = robjects.r['kernel']
  k2 = kernel('modified.daniell', spans//2)
  kernapply2 = robjects.r['kernapply']
  b = flex.double(kernapply2(x, k2, circular=True))

  for e1, e2 in zip(a, b):
    assert approx_equal(e1, e2)
  print "OK"

def test5():
  # compare smoothed pgrams with convolved kernel, even and odd length sequences
  for i in range(2):
    dat = flex.random_double(50+i)

    a = Periodogram(dat, spans=[4,4])
    b = rpgram(dat, spans=robjects.FloatVector(list([4,4])))

    #from matplotlib.pyplot import ion
    #ion()
    #a.plot()
    #b.plot()

    assert approx_equal(a.freq, b.freq)
    assert approx_equal(a.spec, b.spec)

    print "OK"

if __name__=="__main__":
  test1()
  #test2()
  test3()
  test4()
  test5()

  print "OK"
