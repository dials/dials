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
from periodogram import Periodogram
from libtbx.test_utils import approx_equal

class rpgram(object):
  """Calculate a raw periodogram using R"""
  def __init__(self, x, detrend=True):
    dat = robjects.FloatVector(list(x))
    pgram = robjects.r['spec.pgram']
    result = pgram(dat, detrend=detrend, taper=0, fast=False, plot=False)
    self.result = result
    self.spec = flex.double(result.rx2('spec'))
    self.freq = flex.double(result.rx2('freq'))
    return

  def plot(self):
    rplot = robjects.r('plot')
    rplot(self.result, main="Periodogram")


def test1():
  # Compare over a range of lengths from 2 to 1000 with random data
  for i in range(2, 1001):
    dat = flex.random_double(i)
    a = Periodogram(dat)
    b = rpgram(dat)

    assert approx_equal(a.freq, b.freq)
    assert approx_equal(a.spec, b.spec)
    print "OK"

def test2():
  dat = flex.random_double(50)
  a = Periodogram(dat)
  b = rpgram(dat)

  from matplotlib.pyplot import ion
  ion()
  a.plot()
  b.plot()
  from dials.util.command_line import interactive_console; interactive_console(); 1/0 #XXXXX DEBUG

if __name__=="__main__":
  test1()
  #test2()
  print "OK"
