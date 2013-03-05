from __future__ import division
from dials.algorithms.integration import FromBeamVectorToXds

class TestFromBeamVectorToXds(object):
    """Test the FromBeamVectorToXds class"""

    def __init__(self):
        """Initialise the transform"""

        from dials.algorithms.integration import XdsCoordinateSystem
        from math import pi

        self.s0 = ( 0.013141995425357206,
                    0.002199999234194632,
                    1.4504754950989514)
        self.s1 = (-0.01752795848400313,
                   -0.24786554213968193,
                    1.4290948735525306)
        self.m2 = ( 0.999975, -0.001289, -0.006968)
        self.phi = 5.83575672475 * pi / 180

        self.xcs = XdsCoordinateSystem(self.s0, self.s1, self.m2, self.phi)
        self.from_beam_vector_to_xds = FromBeamVectorToXds(
            self.xcs, self.s1, self.phi)

    def test_coordinate_of_s1(self):
        """Ensure that the coordinate of s1 is (0, 0, 0)"""

        eps = 1e-7

        # Get the coordinate at s1
        s_dash = self.s1
        phi_dash = self.phi
        c1, c2, c3 = self.from_beam_vector_to_xds(s_dash, phi_dash)

        # Ensure that it is at the origin
        assert(abs(c1 - 0.0) <= eps)
        assert(abs(c2 - 0.0) <= eps)
        assert(abs(c3 - 0.0) <= eps)

        # Test passed
        print "OK"

    def test_e3_coordinate_approximation(self):

        from scitbx import matrix
        from math import pi
        import random

        eps = 1e-4

        # Select a random rotation from phi
        s_dash = self.s1
        phi_dash = self.phi + (2.0*random.random() - 1.0) * pi / 180

        # Calculate the XDS coordinate, this class uses an approximation
        # for c3 = zeta * (phi' - phi)
        c1, c2, c3 = self.from_beam_vector_to_xds(s_dash, phi_dash)

        # Calculate the true value
        p_star = matrix.col(self.s1) - matrix.col(self.s0)
        scale = 1.0 / p_star.length()
        e3 = matrix.col(self.xcs.e3_axis)
        p_star_rot = p_star.rotate(matrix.col(self.m2), phi_dash - self.phi)
        c3_2 = e3.dot(p_star_rot - p_star) * scale

        # Check the true and approximate value are almost equal to 4dp
        assert(abs(c3 - c3_2) < eps)

        # Test passed
        print "OK"

    def run(self):
        """Run all the tests"""
        self.test_coordinate_of_s1()
        self.test_e3_coordinate_approximation()

if __name__ == '__main__':
    test = TestFromBeamVectorToXds()
    test.run()
