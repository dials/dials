from __future__ import division
import unittest

from dials.geometry import XdsCoordinateSystem

class TestXdsCoordinateSystem(unittest.TestCase):
    """Test the XDS coordinate system class"""

    def setUp(self):
        """Initialise the coordinate system"""

        # Coordinate sensitive to length of vectors, so need to ensure that
        # lengths of both s0 and s1 are equal
        self.s0 = ( 0.013141995425357206, 0.002199999234194632, 1.4504754950989514)
        self.s1 = (-0.01752795848400313, -0.24786554213968193, 1.4290948735525306)
        self.m2 = ( 0.999975, -0.001289, -0.006968)
        self.phi = 5.83575672475
        self.xcs = XdsCoordinateSystem(self.s0, self.s1, self.m2, self.phi)

    def test_axis_length(self):
        """Ensure axes are of unit length"""
        from scitbx import matrix

        # Check all axes are of length 1
        self.assertAlmostEqual(matrix.col(self.xcs.e1_axis).length(), 1.0)
        self.assertAlmostEqual(matrix.col(self.xcs.e2_axis).length(), 1.0)
        self.assertAlmostEqual(matrix.col(self.xcs.e3_axis).length(), 1.0)

    def test_axis_orthogonal(self):
        """Ensure that the following are true:

        e1.s0 = 0, e1.s1 = 0
        e2.s1 = 0, e2.e1 = 0
        e3.e1 = 0, e3.p* = 0

        """
        from scitbx import matrix

        # Get as matrices
        e1 = matrix.col(self.xcs.e1_axis)
        e2 = matrix.col(self.xcs.e2_axis)
        e3 = matrix.col(self.xcs.e3_axis)
        s0 = matrix.col(self.s0)
        s1 = matrix.col(self.s1)

        # Check that e1 is orthogonal to s0 and s1
        self.assertAlmostEqual(e1.dot(s0), 0.0)
        self.assertAlmostEqual(e1.dot(s1), 0.0)

        # Check that e2 is orthogonal to s1 and e1
        self.assertAlmostEqual(e2.dot(s1), 0.0)
        self.assertAlmostEqual(e2.dot(e1), 0.0)

        # Check that e3 is orthogonal to e1 and p* = s1 - s0
        self.assertAlmostEqual(e3.dot(e1), 0.0)
        self.assertAlmostEqual(e3.dot(s1 - s0), 0.0)

if __name__ == '__main__':
    unittest.main()
