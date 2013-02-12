
from detector.coordinate_system import CoordinateSystem

import unittest

class TestDetectorCoordinateSystem(unittest.TestCase):

    def setUp(self):

        from scitbx import matrix

        # Set the parameters needed for the coordinate system
        d1 = matrix.col((0.000000, -0.939693, -0.342020))
        d2 = matrix.col((1.000000,  0.000000,  0.000000))
        d3 = matrix.col((0.000000, -0.342020,  0.939693))
        f = 122.124901
        x0 = 244.836136
        y0 = 320.338531
        px = 0.172000
        py = 0.172000

        # Create the coordinate system
        self.cs = CoordinateSystem(d1, d2, d3, f, (x0, y0), (px, py))

    def test_origin_pixel_coordinates(self):

        # Set the wavelength
        wavelength = 0.689400

        # Set the beam vector to directed towards the detector origin
        s = (1.0 / wavelength) * (self.cs.f * self.cs.d3).normalize()

        # Get the x, y pixel coordinates
        x, y = self.cs.from_beam_vector(s)

        # Check that they are equal to the detector origin
        self.assertAlmostEqual(x, self.cs.x0)
        self.assertAlmostEqual(y, self.cs.y0)

    def test_away_from_detector_pixel_coordinates(self):

        # Set the wavelength
        wavelength = 0.689400

        # Set the beam vector to directed away from the detector origin
        s = - (1.0 / wavelength) * (self.cs.f * self.cs.d3).normalize()

        # Get the x, y pixel coordinates
        x, y = self.cs.from_diffracted_beam_vector(s)

        # Check that they are equal to None
        self.assertAlmostEqual(x, None)
        self.assertAlmostEqual(y, None)

    def test_origin_laboratory_coordinates(self):

        # Get the origin in laboratory coordinates
        origin = self.cs.to_laboratory_coordinate(self.cs.x0, self.cs.y0)

        # Calculate the actual origin
        real_origin = self.cs.f * self.cs.d3

        # Check they match
        self.assertAlmostEqual(origin, real_origin)
