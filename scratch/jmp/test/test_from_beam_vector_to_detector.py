import unittest

from dials.geometry.transform import FromBeamVectorToDetector

class TestFromBeamVectorToDetector(unittest.TestCase):
    """Test the beam vector to detector transform class"""

    def setUp(self):
        """Initialise the transform"""

        from dials.geometry import DetectorCoordinateSystem

        # The input parameters
        self.x_axis = (0.000000, -0.939693, -0.342020)
        self.y_axis = (1.000000,  0.000000,  0.000000)
        self.normal = (0.000000, -0.342020,  0.939693)
        self.pixel_size = (0.172000, 0.172000)
        self.origin = (244.836136, 320.338531)
        self.distance = 122.124901

        # Create the transform object
        self.from_beam_vector_to_detector = FromBeamVectorToDetector(
            DetectorCoordinateSystem(
                self.x_axis,
                self.y_axis,
                self.normal),
            self.pixel_size,
            self.origin,
            self.distance)

    def test_origin_pixel_coordinates(self):
        """Ensure that the detector coordinates of the beam vector at the
        detector origin equals the detector origin."""

        from scitbx import matrix

        # Set the wavelength
        wavelength = 0.689400

        # Set the beam vector to directed towards the detector origin
        s = (1.0 / wavelength) * (self.distance *
            matrix.col(self.normal)).normalize()

        # Get the x, y pixel coordinates
        x, y = self.from_beam_vector_to_detector.apply(s)

        # Check that they are equal to the detector origin
        self.assertAlmostEqual(x, self.origin[0])
        self.assertAlmostEqual(y, self.origin[1])

    def test_away_from_detector_pixel_coordinates(self):
        """For a beam vector pointing away from the detector, ensure that the
        transform raises an exception indicating the detector coordinate of the
        beam vector cannot be calculated."""

        from scitbx import matrix

        # Set the wavelength
        wavelength = 0.689400

        # Set the beam vector to directed away from the detector origin
        s = - (1.0 / wavelength) * (self.distance *
            matrix.col(self.normal)).normalize()

        # Check that they are equal to None
        with self.assertRaises(RuntimeError):
            self.from_beam_vector_to_detector.apply(s)

if __name__ == '__main__':
    unittest.main()
