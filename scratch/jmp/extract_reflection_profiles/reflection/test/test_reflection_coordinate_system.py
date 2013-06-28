from __future__ import division

from reflection.coordinate_system import CoordinateSystem

import unittest

class TestReflectionCoordinateSystem(unittest.TestCase):

    def setUp(self):

        from scitbx import matrix

        wavelength = 0.689400
        s0 = matrix.col(( 0.013142,  0.002200,  1.450476))
        s1 = matrix.col((-0.017520, -0.247753,  1.428447))
        m2 = matrix.col(( 0.999975, -0.001289, -0.006968))
        phi = 5.835757

        s0 = s0.normalize() / wavelength
        s1 = s1.normalize() / wavelength

        self.cs = CoordinateSystem(s0, s1, m2, phi)

    def test_axis_length(self):

        # Check that all axes are of length 1
        self.assertAlmostEqual(self.cs.e1.length(), 1.0)
        self.assertAlmostEqual(self.cs.e2.length(), 1.0)
        self.assertAlmostEqual(self.cs.e3.length(), 1.0)

    def test_axis_orthogonal(self):

        # Check that e1 is orthogonal to s0 and s1
        self.assertAlmostEqual(self.cs.e1.dot(self.cs.s0), 0.0)
        self.assertAlmostEqual(self.cs.e1.dot(self.cs.s1), 0.0)

        # Check that e2 is orthogonal to s1 and e1
        self.assertAlmostEqual(self.cs.e2.dot(self.cs.s1), 0.0)
        self.assertAlmostEqual(self.cs.e2.dot(self.cs.e1), 0.0)

        # Check that e3 is orthogonal to e1 and p*
        self.assertAlmostEqual(self.cs.e3.dot(self.cs.e1), 0.0)
        self.assertAlmostEqual(self.cs.e3.dot(self.cs.s1 - self.cs.s0), 0.0)

    def test_coordinate_of_s1(self):

        # Get the coordinate at s1
        s_dash = self.cs.s1
        phi_dash = self.cs.phi
        x1, x2, x3 = self.cs.calculate_coordinate(s_dash, phi_dash)

        # Ensure that it is at the origin
        self.assertAlmostEqual(x1, 0.0)
        self.assertAlmostEqual(x2, 0.0)
        self.assertAlmostEqual(x3, 0.0)

if __name__ == '__main__':
    unittest.main()
