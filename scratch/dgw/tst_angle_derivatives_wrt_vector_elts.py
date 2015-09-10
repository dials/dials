#!/usr/bin/env dials.python

"""Test analytical expression for the partial derivatives of an angle
between two vectors with respect to each element of the vectors"""

from __future__ import division
from scitbx import matrix
from random import uniform
from math import sin, cos, acos
from libtbx.test_utils import approx_equal

def dgamma_da1(a, b):
  ab = a.length()*b.length()
  gamma = acos(a.dot(b) / (ab))
  a1 = a.elems[0]
  b1 = b.elems[0]

  return -1. * (b1/ab - a1*cos(gamma)/a.length()**2) / sin(gamma)

def dgamma_da1_fd(a, b, delta=1.e-7):

  half_delta = matrix.col((delta/2., 0, 0))
  a_fwd = a + half_delta
  a_rev = a - half_delta

  gamma_fwd = acos(a_fwd.dot(b) / (a_fwd.length()*b.length()))
  gamma_rev = acos(a_rev.dot(b) / (a_rev.length()*b.length()))

  return (gamma_fwd - gamma_rev) / delta

def test():
  # two random vectors
  vec_a = matrix.col((uniform(-20,20), uniform(-20,20), uniform(-20,20)))
  vec_b = matrix.col((uniform(-20,20), uniform(-20,20), uniform(-20,20)))

  a = vec_a.length()
  b = vec_b.length()

  a_dot_b = vec_a.dot(vec_b)

  # angle
  gamma = acos(a_dot_b / (a*b))

  dg_da1 = dgamma_da1(vec_a, vec_b)

  dg_da1_fd = dgamma_da1_fd(vec_a, vec_b)

  assert approx_equal(dg_da1, dg_da1_fd)

if __name__ == '__main__':

  for i in range(100): test()




