
from __future__ import division
from numpy.random import poisson

def run(mean, n1, n2):

  D = []
  for i in range(n2):
    print i
    y = list(poisson(mean, n1))
    m = sum(y) / len(y)
    v = sum(yy*yy for yy in y) / len(y) - m*m
    assert(m >= 0)
    assert(v >= 0)
    if m > 0:
      D.append(v / m)

  print len(D)
  from matplotlib import pylab
  pylab.hist(D, bins=50)
  pylab.show()

run(0.01, 10000, 100000)
