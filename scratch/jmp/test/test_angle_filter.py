import unittest

from dials.spot_prediction import angle_filter

class TestAngleFilter(unittest.TestCase):
    """Test the dials.spot_prediction.angle_filter function"""
    
    def setUp(self):
        pass
        
    def test_radians(self):
        """Test the angle filter with ranges between -pi and pi"""
        import random
        from math import pi
        random_angle = lambda: random.random() * 2 * pi - pi

        angle_range = (random_angle(), random_angle())
        angle_range = (min(angle_range), max(angle_range))
        
        for i in range(0, 1000):
            angle = random_angle()
            within = angle_filter(angle, angle_range)
            if (angle_range[0] <= angle <= angle_range[1]):
                self.assertTrue(within)
            else:
                self.assertFalse(within)

    def test_degrees(self):
        """Test the angle filter with ranges between -180 and 180"""        
        import random
        from math import pi
        random_angle = lambda: random.random() * 360 - 180

        angle_range = (random_angle(), random_angle())
        angle_range = (min(angle_range), max(angle_range))
        
        for i in range(0, 1000):
            angle = random_angle()
            within = angle_filter(angle, angle_range, True)
            if (angle_range[0] <= angle <= angle_range[1]):
                self.assertTrue(within)
            else:
                self.assertFalse(within)

if __name__ == '__main__':
    unittest.main()
