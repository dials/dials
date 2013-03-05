from __future__ import division
from dials.algorithms.integration import FromXdsE3ToPhi

class TestFromXdsE3ToPhi(object):
    """Test the FromXdsE3ToPhi class"""

    def __init__(self):
        """Initialise the transform"""

        from dials.algorithms.integration import XdsCoordinateSystem
        from dials.algorithms.integration import FromBeamVectorToXds
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
        self.from_xds_e3_to_phi = FromXdsE3ToPhi(self.xcs.zeta, self.phi)
        self.from_beam_vector_to_xds = FromBeamVectorToXds(self.xcs, self.s1,
                                                           self.phi)

    def test_forward_and_backward(self):

        from scitbx import matrix
        from math import pi
        import random
        eps = 1e-7

        # Set the parameters
        s_dash = self.s1
        min_shift = -5.0 * pi / 180.0
        max_shift = +5.0 * pi / 180.0
        mod_2pi_shift = lambda: int(4 * random.random() - 2)
        range_shift = max_shift - min_shift
        random_shift = lambda: min_shift + random.random() * range_shift

        # Loop a number of times
        num = 1000
        for i in range(num):

            # Create a rotation angle
            phi_dash = self.phi + mod_2pi_shift() * 2*pi + random_shift()

            # Calculate the XDS coordinate of the vector
            c1, c2, c3 = self.from_beam_vector_to_xds(s_dash, phi_dash)

            # Calculate the beam vector from the XDS coordinate
            phi_dash_2 = self.from_xds_e3_to_phi(c3)

            # Check the vectors are almost equal
            assert(abs(phi_dash_2 - phi_dash) <= eps)

        # Test passed
        print "OK"

    def run(self):
        """Run all the tests"""
        self.test_forward_and_backward()

if __name__ == '__main__':
    test = TestFromXdsE3ToPhi()
    test.run()
