import unittest

from dials.geometry.transform import FromDetectorToBeamVector

class TestFromDetectorToBeamVector(unittest.TestCase):
    """Test the dials.geometry.transform.FromDetectorToBeamVector class"""

    def setUp(self):
        """Initialise the transform"""

        from dials.geometry import DetectorCoordinateSystem
        from dials.geometry.transform import FromBeamVectorToDetector

        self.x_axis = (0.000000, -0.939693, -0.342020)
        self.y_axis = (1.000000,  0.000000,  0.000000)
        self.normal = (0.000000, -0.342020,  0.939693)
        self.pixel_size = (0.172000, 0.172000)
        self.origin = (244.836136, 320.338531)
        self.distance = 122.124901
        self.sxy = (117.588714455, 311.621428845)

        self.from_detector_to_beam_vector = FromDetectorToBeamVector(
            DetectorCoordinateSystem(
                self.x_axis,
                self.y_axis,
                self.normal),
            self.pixel_size,
            self.origin, self.distance)

        self.from_beam_vector_to_detector = FromBeamVectorToDetector(
            DetectorCoordinateSystem(
                self.x_axis,
                self.y_axis,
                self.normal),
            self.pixel_size,
            self.origin, self.distance)

    def test_detector_origin(self):
        """Test the coordinate of the detector origin"""

        from scitbx import matrix

        # Get the origin in laboratory coordinates
        origin = matrix.col(self.from_detector_to_beam_vector.apply(self.origin))

        # Calculate the actual origin
        real_origin = self.distance * matrix.col(self.normal).normalize()

        # Check they match
        self.assertAlmostEqual(origin, real_origin)

    def test_forward_and_reverse_transform(self):
        """Test the forward and reverse detector -> beam vector transforms"""

        from scitbx import matrix
        import random

        # Set the shift parameters
        min_shift = -5
        max_shift = +5
        range_shift = max_shift - min_shift
        random_shift = lambda: min_shift + random.random() * range_shift

        # Loop a number of times
        num = 1000
        for i in range(num):

            # Create a beam vector
            xy = matrix.col(self.sxy) + matrix.col((random_shift(),
                                                    random_shift()))

            # Calculate the XDS coordinate of the vector
            s = self.from_detector_to_beam_vector.apply(xy)

            # Calculate the beam vector from the XDS coordinate
            xy_2 = self.from_beam_vector_to_detector.apply(s)

            # Check the vectors are almost equal
            self.assertAlmostEqual(xy, matrix.col(xy_2))

if __name__ == '__main__':
    unittest.main()
