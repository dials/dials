from __future__ import division
import unittest

from dials.geometry.transform import FromBeamVectorToXds

class TestFromBeamVectorToXds(unittest.TestCase):
    """Test the dials.geometry.transform.FromBeamVectorToXds class"""

    def setUp(self):
        """Initialise the transform"""

        from dials.geometry import XdsCoordinateSystem

        self.s0 = ( 0.013141995425357206, 0.002199999234194632, 1.4504754950989514)
        self.s1 = (-0.01752795848400313, -0.24786554213968193, 1.4290948735525306)
        self.m2 = ( 0.999975, -0.001289, -0.006968)
        self.phi = 5.83575672475
        self.xcs = XdsCoordinateSystem(self.s0, self.s1, self.m2, self.phi)
        self.from_beam_vector_to_xds = FromBeamVectorToXds(
            self.xcs, self.s1, self.phi)

    def test_coordinate_of_s1(self):
        """Ensure that the coordinate of s1 is (0, 0, 0)"""

        # Get the coordinate at s1
        s_dash = self.s1
        phi_dash = self.phi
        c1, c2, c3 = self.from_beam_vector_to_xds.apply(s_dash, phi_dash)

        # Ensure that it is at the origin
        self.assertAlmostEqual(c1, 0.0)
        self.assertAlmostEqual(c2, 0.0)
        self.assertAlmostEqual(c3, 0.0)

    def test_e3_coordinate_approximation(self):

        from scitbx import matrix
        from math import pi
        import random

        # Select a random rotation from phi
        s_dash = self.s1
        phi_dash = self.phi + 2.0*random.random() - 1.0

        # Calculate the XDS coordinate, this class uses an approximation
        # for c3 = zeta * (phi' - phi)
        c1, c2, c3 = self.from_beam_vector_to_xds.apply(s_dash, phi_dash)

        # Calculate the true value
        p_star = matrix.col(self.s1) - matrix.col(self.s0)
        scale = 180.0 / (p_star.length() * pi)
        e3 = matrix.col(self.xcs.e3_axis)
        p_star_rot = p_star.rotate(matrix.col(self.m2), phi_dash - self.phi, deg=True)
        c3_2 = e3.dot(p_star_rot - p_star) * scale

        # Check the true and approximate value are almost equal to 4dp
        self.assertAlmostEqual(c3, c3_2, places=4)

if __name__ == '__main__':
    unittest.main()
