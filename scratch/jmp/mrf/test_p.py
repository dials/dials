
from __future__ import division

def P(mu, x):
  from math import exp, factorial
  return exp(-mu)*(mu**x)/factorial(x)

def prd(x):
  r = 1
  for xx in x:
    r *= xx
  return r

if __name__ == '__main__':
  from numpy.random import poisson

  x = list(poisson(1, 50))

  p = [P(1,xx) for xx in x]

  print p
