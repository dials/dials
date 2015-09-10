#!/usr/bin/env dials.python

"""Test analytical expression for the partial derivatives of an angle
between two vectors with respect to each element of the vectors"""

from __future__ import division
from scitbx import matrix
from random import uniform
from math import sin, cos, acos
from libtbx.test_utils import approx_equal

from dials.algorithms.refinement.refinement_helpers import \
  AngleDerivativeWrtVectorElts

class FDAngleDerivativeWrtVectorElts(object):
  '''Given two vectors, a and b, calculate the derivative of the angle gamma
  between them with respect to any of the elements a_1, a_2, a_3, b_1, b_2
  and b_3 using a finite difference approximation'''

  def __init__(self, a, b, delta=1.e-7):

    self._vec_a = a
    self._vec_b = b
    self._a = a.length()
    self._b = b.length()
    self._delta = delta

    return

  def dgamma_da_elt(self, i):
    '''Return the derivative of gamma with respect to the ith element of a'''

    half_delta_shift = [0., 0., 0.]
    half_delta_shift[i] = self._delta/2.
    half_delta_shift = matrix.col(half_delta_shift)

    a_fwd = self._vec_a + half_delta_shift
    a_rev = self._vec_a - half_delta_shift

    gamma_fwd = acos(a_fwd.dot(self._vec_b) / (a_fwd.length()*self._b))
    gamma_rev = acos(a_rev.dot(self._vec_b) / (a_rev.length()*self._b))

    return (gamma_fwd - gamma_rev) / self._delta

  def dgamma_db_elt(self, i):
    '''Return the derivative of gamma with respect to the ith element of b'''

    half_delta_shift = [0., 0., 0.]
    half_delta_shift[i] = self._delta/2.
    half_delta_shift = matrix.col(half_delta_shift)

    b_fwd = self._vec_b + half_delta_shift
    b_rev = self._vec_b - half_delta_shift

    gamma_fwd = acos(self._vec_a.dot(b_fwd) / (b_fwd.length()*self._a))
    gamma_rev = acos(self._vec_a.dot(b_rev) / (b_rev.length()*self._a))

    return (gamma_fwd - gamma_rev) / self._delta

  def dgamma_da_1(self): return self.dgamma_da_elt(0)

  def dgamma_da_2(self): return self.dgamma_da_elt(1)

  def dgamma_da_3(self): return self.dgamma_da_elt(2)

  def dgamma_db_1(self): return self.dgamma_db_elt(0)

  def dgamma_db_2(self): return self.dgamma_db_elt(1)

  def dgamma_db_3(self): return self.dgamma_db_elt(2)

def test():

  # Two random vectors
  vec_a = matrix.col((uniform(-20,20), uniform(-20,20), uniform(-20,20)))
  vec_b = matrix.col((uniform(-20,20), uniform(-20,20), uniform(-20,20)))

  # The test may fail if the angle between the vectors is small or close to pi
  # (see gradient of acos(x) for x close to 1 or -1) but we are not interested
  # in such cases anyway. Skip test if the acute angle is less than 1 degree
  if vec_a.accute_angle(vec_b, deg=True) < 1: return False

  a = vec_a.length()
  b = vec_b.length()
  a_dot_b = vec_a.dot(vec_b)
  gamma = acos(a_dot_b / (a*b))

  # Analytical
  dg = AngleDerivativeWrtVectorElts(vec_a, vec_b)

  # FD
  dgFD = FDAngleDerivativeWrtVectorElts(vec_a, vec_b)
  assert approx_equal(dg.dgamma_da_1(), dgFD.dgamma_da_1())
  assert approx_equal(dg.dgamma_da_2(), dgFD.dgamma_da_2())
  assert approx_equal(dg.dgamma_da_3(), dgFD.dgamma_da_3())
  assert approx_equal(dg.dgamma_db_1(), dgFD.dgamma_db_1())
  assert approx_equal(dg.dgamma_db_2(), dgFD.dgamma_db_2())
  assert approx_equal(dg.dgamma_db_3(), dgFD.dgamma_db_3())

  return True

if __name__ == '__main__':

  results = [test() for i in xrange(10000)]

  print "{0} successes".format(results.count(True))
  print "{0} skipped tests".format(results.count(False))
