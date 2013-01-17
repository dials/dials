import unittest

from dials.geometry.transform import FromDetectorToXds
from dials.geometry.transform import FromXdsToDetector

class TestFromDetectorToXds(unittest.TestCase):
    """Test the dials.geometry.transform FromDetectorToXds class"""

    def setUp(self):
        """Initialise the transform"""
        
        from dials.geometry.transform import FromBeamVectorToDetector
        from dials.geometry.transform import FromBeamVectorToXds
        from dials.geometry.transform import FromDetectorToBeamVector
        from dials.geometry.transform import FromXdsToBeamVector
        from dials.geometry import XdsCoordinateSystem
        from dials.geometry import DetectorCoordinateSystem
        from scitbx import matrix

        # Set the parameters
        s0 = ( 0.013141995425357206, 0.002199999234194632, 1.4504754950989514)
        s1 = (-0.01752795848400313, -0.24786554213968193, 1.4290948735525306)
        m2 = ( 0.999975, -0.001289, -0.006968)
        phi = 5.83575672475
        x_axis = (0.000000, -0.939693, -0.342020)
        y_axis = (1.000000,  0.000000,  0.000000)
        normal = (0.000000, -0.342020,  0.939693)
        pixel_size = (0.172000, 0.172000)
        origin = (244.836136, 320.338531)
        distance = 122.124901
        sxy = (117.588714455, 311.621428845)
        wavelength = 1.0 / matrix.col(s1).length()
        
        # Create the component transforms
        dcs = DetectorCoordinateSystem(x_axis, y_axis, normal)
        xcs = XdsCoordinateSystem(s0, s1, m2, phi)
        from_beam_vector_to_detector = FromBeamVectorToDetector(
            dcs, pixel_size, origin, distance)
        from_detector_to_beam_vector = FromDetectorToBeamVector(
            dcs, pixel_size, origin, distance)
        from_beam_vector_to_xds = FromBeamVectorToXds(xcs, s1, phi)
        from_xds_to_beam_vector = FromXdsToBeamVector(xcs, s1)

        # Create the forward and reverse transforms
        self.from_detector_to_xds = FromDetectorToXds(
                                        from_detector_to_beam_vector,
                                        from_beam_vector_to_xds,
                                        wavelength)
        self.from_xds_to_detector = FromXdsToDetector(
                                        from_xds_to_beam_vector,
                                        from_beam_vector_to_detector)

        self.sxy = sxy
        self.s1 = s1
        self.phi = phi

    def test_forward_and_reverse_transform(self):
        """Test the forward and reverse Detector -> XDS transforms Create
        a detector coordinare, transform it to XDS and then transform back. 
        The new value should be equal to the original value."""
        
        from scitbx import matrix
        import random

        # Set the shift parameters
        phi_dash = self.phi
        min_shift = -0.5
        max_shift = +0.5
        range_shift = max_shift - min_shift
        random_shift = lambda: min_shift + random.random() * range_shift

        # Loop a number of times
        num = 1000
        for i in range(num):
                
            # Create a beam vector
            xy = matrix.col(self.sxy) + matrix.col((random_shift(),
                                                    random_shift()))
        
            # Calculate the XDS coordinate of the vector
            c1, c2, c3 = self.from_detector_to_xds.apply(xy, phi_dash)

            # Calculate the beam vector from the XDS coordinate
            xy_2 = self.from_xds_to_detector.apply((c1, c2, c3))

            # Check the vectors are almost equal
            self.assertAlmostEqual(xy, matrix.col(xy_2))

        
if __name__ == '__main__':
    unittest.main()
