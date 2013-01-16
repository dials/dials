import unittest

from dials.geometry import ReciprocalLatticeCoordinateSystem

class TestReciprocalLatticeCoordinateSystem(unittest.TestCase):
    """Test the dials.geometry.ReciprocalLatticeCoordinateSystem class"""
    
    def setUp(self):
        pass
    
    def test_data(self):
        """Test the member data"""
        
        from scitbx import matrix
        
        b1 = (33.915398, -36.200676, -14.127665)
        b2 = (38.855934,  31.328993,  13.001726)
        b3 = (-0.544143, -19.192179,  47.871685)
        ub = matrix.sqr(b1 + b2 + b3).inverse()
        
        b1_star = (ub[0], ub[3], ub[6])
        b2_star = (ub[1], ub[4], ub[7])
        b3_star = (ub[2], ub[5], ub[8])
        rcs = ReciprocalLatticeCoordinateSystem(b1_star, b2_star, b3_star)
        
        self.assertEqual(rcs.b1_star_axis, b1_star)
        self.assertEqual(rcs.b2_star_axis, b2_star)
        self.assertEqual(rcs.b3_star_axis, b3_star)
     
    def test_init_from_ub(self):
        """Test the initialisation from UB matrix"""
    
        from scitbx import matrix
        
        b1 = (33.915398, -36.200676, -14.127665)
        b2 = (38.855934,  31.328993,  13.001726)
        b3 = (-0.544143, -19.192179,  47.871685)
        ub = matrix.sqr(b1 + b2 + b3).inverse()
        
        b1_star = (ub[0], ub[3], ub[6])
        b2_star = (ub[1], ub[4], ub[7])
        b3_star = (ub[2], ub[5], ub[8])
        rcs = ReciprocalLatticeCoordinateSystem(ub)
        
        self.assertEqual(rcs.b1_star_axis, b1_star)
        self.assertEqual(rcs.b2_star_axis, b2_star)
        self.assertEqual(rcs.b3_star_axis, b3_star)
    
    def test_ub_from_rcs(self):
        """Test the conversion to UB matrix method"""
        
        from scitbx import matrix
        
        b1 = (33.915398, -36.200676, -14.127665)
        b2 = (38.855934,  31.328993,  13.001726)
        b3 = (-0.544143, -19.192179,  47.871685)
        ub = matrix.sqr(b1 + b2 + b3).inverse()
        rcs = ReciprocalLatticeCoordinateSystem(ub)
        ub2 = rcs.to_ub_matrix()
        self.assertEqual(tuple(ub), ub2)
        
if __name__ == '__main__':
    unittest.main()
