import unittest

from dials.equipment import Goniometer

class TestGoniometer(unittest.TestCase):
    """Test the dials.equipment.Goniometer class."""
    
    def setUp(self):
        """Initialise the goniometer class."""
        
        self.rotation_axis = (0.999975, -0.001289, -0.006968)
        self.starting_angle = 1.0
        self.oscillation_range = 1.0
        self.starting_frame = 1
        
        self.goniometer = Goniometer(self.rotation_axis,
                                     self.starting_angle,
                                     self.oscillation_range,
                                     self.starting_frame)

    
    def test_data(self):
        """Test the class member data."""
    
        self.assertEqual(self.goniometer.rotation_axis, self.rotation_axis)                                
        self.assertEqual(self.goniometer.starting_angle, self.starting_angle)
        self.assertEqual(self.goniometer.oscillation_range, self.oscillation_range)
        self.assertEqual(self.goniometer.starting_frame, self.starting_frame)
        
    def test_frame_and_angle_conversion(self):
        """Test the methods for calculating the angle from the frame number and
        the frame number from the angle. Ensure that giving one function the
        output of the other gives the initial value."""
        
        import random
        n_frames = 1000
        min_frame = -1000
        max_frame = 1000
        random_frame = lambda: min_frame + random.random() * (max_frame - min_frame)
        frames = [random_frame() for x in range(n_frames)]
        for f in frames:
            phi = self.goniometer.get_angle_from_frame(f)
            f2  = self.goniometer.get_frame_from_angle(phi)
            self.assertAlmostEqual(f, f2)
        
    
        
if __name__ == '__main__':
    unittest.main()
