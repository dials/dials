"""Test analytical expression for the partial derivatives of an angle
between two vectors with respect to each element of the vectors"""

from __future__ import absolute_import, division, print_function

import math
import random

from libtbx.test_utils import approx_equal
from scitbx import matrix
from scitbx.math import angle_derivative_wrt_vectors


class FDAngleDerivativeWrtVectorElts(object):
    """Given two vectors, u and v, calculate the derivative of the angle theta
    between them with respect to any of the elements u_1, u_2, u_3, v_1, v_2
    and v_3 using a finite difference approximation"""

    def __init__(self, u, v, delta=1.0e-7):
        self._vec_u = u
        self._vec_v = v
        self._u = u.length()
        self._v = v.length()
        self._delta = delta

    def dtheta_du_elt(self, i):
        """Return the derivative of theta with respect to the ith element of a"""

        half_delta_shift = [0.0, 0.0, 0.0]
        half_delta_shift[i] = self._delta / 2.0
        half_delta_shift = matrix.col(half_delta_shift)

        u_fwd = self._vec_u + half_delta_shift
        u_rev = self._vec_u - half_delta_shift

        theta_fwd = math.acos(u_fwd.dot(self._vec_v) / (u_fwd.length() * self._v))
        theta_rev = math.acos(u_rev.dot(self._vec_v) / (u_rev.length() * self._v))

        return (theta_fwd - theta_rev) / self._delta

    def dtheta_dv_elt(self, i):
        """Return the derivative of theta with respect to the ith element of b"""

        half_delta_shift = [0.0, 0.0, 0.0]
        half_delta_shift[i] = self._delta / 2.0
        half_delta_shift = matrix.col(half_delta_shift)

        v_fwd = self._vec_v + half_delta_shift
        v_rev = self._vec_v - half_delta_shift

        theta_fwd = math.acos(self._vec_u.dot(v_fwd) / (v_fwd.length() * self._u))
        theta_rev = math.acos(self._vec_u.dot(v_rev) / (v_rev.length() * self._u))

        return (theta_fwd - theta_rev) / self._delta

    def dtheta_du_1(self):
        return self.dtheta_du_elt(0)

    def dtheta_du_2(self):
        return self.dtheta_du_elt(1)

    def dtheta_du_3(self):
        return self.dtheta_du_elt(2)

    def dtheta_dv_1(self):
        return self.dtheta_dv_elt(0)

    def dtheta_dv_2(self):
        return self.dtheta_dv_elt(1)

    def dtheta_dv_3(self):
        return self.dtheta_dv_elt(2)


def _test():
    # Two random vectors
    vec_a = matrix.col(
        (random.uniform(-20, 20), random.uniform(-20, 20), random.uniform(-20, 20))
    )
    vec_b = matrix.col(
        (random.uniform(-20, 20), random.uniform(-20, 20), random.uniform(-20, 20))
    )

    # The test may fail if the angle between the vectors is small or close to pi
    # (see gradient of acos(x) for x close to 1 or -1) but we are not interested
    # in such cases anyway. Skip test if the acute angle is less than 1 degree
    if vec_a.accute_angle(vec_b, deg=True) < 1:
        return False

    # Analytical
    def dangle(u, v):
        return [matrix.col(e) for e in angle_derivative_wrt_vectors(u, v)]

    dgamma_da, dgamma_db = dangle(vec_a, vec_b)

    # FD
    dgFD = FDAngleDerivativeWrtVectorElts(vec_a, vec_b)
    assert approx_equal(
        dgamma_da,
        matrix.col([dgFD.dtheta_du_1(), dgFD.dtheta_du_2(), dgFD.dtheta_du_3()]),
    )
    assert approx_equal(
        dgamma_db,
        matrix.col([dgFD.dtheta_dv_1(), dgFD.dtheta_dv_2(), dgFD.dtheta_dv_3()]),
    )

    return True


def test_should_succeed_in_95_percent_of_cases():
    ntests = 1000
    results = [_test() for i in range(ntests)]
    nsuccess = results.count(True)
    assert nsuccess > 0.95 * ntests
