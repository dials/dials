from __future__ import division
import unittest

from dials.geometry.transform import FromXdsToBeamVector

class TestFromXdsToBeamVector(unittest.TestCase):
    """Test the dials.geometry.transform.FromXdsToBeamVector class"""

    def setUp(self):
        """Initialise the transform"""

        from dials.geometry.transform import FromBeamVectorToXds
        from dials.geometry import XdsCoordinateSystem

        self.s0 = ( 0.013141995425357206, 0.002199999234194632, 1.4504754950989514)
        self.s1 = (-0.01752795848400313, -0.24786554213968193, 1.4290948735525306)
        self.m2 = ( 0.999975, -0.001289, -0.006968)
        self.phi = 5.83575672475

        self.xcs = XdsCoordinateSystem(self.s0, self.s1, self.m2, self.phi)
        self.from_beam_vector_to_xds = FromBeamVectorToXds(self.xcs, self.s1, self.phi)
        self.from_xds_to_beam_vector = FromXdsToBeamVector(self.xcs, self.s1)

    def test_xds_origin(self):
        """Test the beam vector at the XDS origin is equal to s1."""

        s_dash = self.from_xds_to_beam_vector.apply((0, 0, 0))
        self.assertAlmostEqual(s_dash, self.s1)

    def test_far_out_coordinates(self):
        """Test some large coordinates, 1 valid and the other invalid. (i.e.
        a coordinate that cannot be mapped onto the ewald sphere)."""

        from scitbx import matrix
        from math import pi
        wavelength = 1.0 / matrix.col(self.s1).length()
        scale = 180.0 / (matrix.col(self.s1).length() * pi)

        # Settting c2 and c3 to zero
        c2 = 0
        c3 = 0

        # A large value which is still valid
        c1 = (matrix.col(self.s1).length() - 0.00001) * scale
        s_dash = self.from_xds_to_beam_vector.apply((c1, c2, c3))

        # A large value which is raises an exception
        with self.assertRaises(RuntimeError):
            c1 = (matrix.col(self.s1).length() + 0.00001) * scale
            s_dash = self.from_xds_to_beam_vector.apply((c1, c2, c3))

    def test_forward_and_reverse_transform(self):
        """Test the forward and reverse Beam Vector -> XDS transforms Create
        a beam vector, transform it to XDS and then transform back. The new
        value should be equal to the original value."""

        from scitbx import matrix
        import random

        # Set the parameters
        phi_dash = self.phi
        min_shift = -0.5
        max_shift = +0.5
        range_shift = max_shift - min_shift
        random_shift = lambda: min_shift + random.random() * range_shift

        # Loop a number of times
        num = 1000
        for i in range(num):

            # Create a beam vector
            s_dash = matrix.col(self.s1) + matrix.col((random_shift(),
                                                       random_shift(),
                                                       random_shift()))
            s_dash = s_dash.normalize() * matrix.col(self.s1).length()

            # Calculate the XDS coordinate of the vector
            c1, c2, c3 = self.from_beam_vector_to_xds.apply(s_dash, phi_dash)

            # Calculate the beam vector from the XDS coordinate
            s_dash_2 = self.from_xds_to_beam_vector.apply((c1, c2, c3))

            # Check the vectors are almost equal
            self.assertAlmostEqual(s_dash, matrix.col(s_dash_2))


if __name__ == '__main__':
    unittest.main()
