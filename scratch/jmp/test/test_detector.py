import unittest

from dials.equipment import Detector

class TestDetector(unittest.TestCase):
    """Test the dials.equipment.Detector class"""
    def setUp(self):
        pass
    
    def test_data(self):
        """Test the member data"""    
        x_axis = (0.000000, -0.939693, -0.342020)
        y_axis = (1.000000,  0.000000,  0.000000)
        normal = (0.000000, -0.342020,  0.939693)
        origin = (244.836136, 320.338531)
        pixel_size = (0.172000, 0.172000)
        size = (487, 619)
        distance = 122.124901
    
        detector = Detector(x_axis, y_axis, normal, 
                            origin, pixel_size, size, distance)
                            
        self.assertEqual(detector.x_axis, x_axis)
        self.assertEqual(detector.y_axis, y_axis)
        self.assertEqual(detector.normal, normal)       
        self.assertEqual(detector.origin, origin)
        self.assertEqual(detector.pixel_size, pixel_size)
        self.assertEqual(detector.size, size)
        self.assertEqual(detector.distance, distance)
                
if __name__ == '__main__':
    unittest.main()
