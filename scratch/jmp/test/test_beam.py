from __future__ import division
import unittest

from dials.equipment import Beam

class TestBeam(unittest.TestCase):
    """Test the dials.equipment.Beam class"""

    def setUp(self):
        pass

    def test_data(self):
        """Test the member data"""
        direction = (0.013142, 0.002200, 1.450476)
        wavelength = 0.689400
        beam = Beam(direction, wavelength)
        self.assertEqual(beam.direction, direction)
        self.assertEqual(beam.wavelength, wavelength)

if __name__ == '__main__':
    unittest.main()
