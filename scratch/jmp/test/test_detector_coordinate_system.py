import unittest

from dials.geometry import DetectorCoordinateSystem

class TestDetectorCoordinateSystem(unittest.TestCase):
    """Test the dials.geometry.DetectorCoordinateSystem class"""

    def setUp(self):
        pass

    def test_data(self):
        """Test the member data"""

        x_axis = (0.000000, -0.939693, -0.342020)
        y_axis = (1.000000,  0.000000,  0.000000)
        normal = (0.000000, -0.342020,  0.939693)
        dcs = DetectorCoordinateSystem(x_axis,
                                       y_axis,
                                       normal)

        self.assertEqual(dcs.x_axis, x_axis)
        self.assertEqual(dcs.y_axis, y_axis)
        self.assertEqual(dcs.normal, normal)

    def test_normal(self):
        """Check that when only given the x and y axes, the calculated normal
        is orthogonal to both x and y axes and is normalized"""

        from scitbx import matrix
        x_axis = (0.000000, -0.939693, -0.342020)
        y_axis = (1.000000,  0.000000,  0.000000)
        dcs = DetectorCoordinateSystem(x_axis,
                                       y_axis)

        normal = matrix.col(dcs.normal)
        x_axis = matrix.col(dcs.x_axis)
        y_axis = matrix.col(dcs.y_axis)
        self.assertAlmostEqual(normal.dot(x_axis), 0.0)
        self.assertAlmostEqual(normal.dot(y_axis), 0.0)
        self.assertAlmostEqual(normal.length(), 1.0)


if __name__ == '__main__':
    unittest.main()
