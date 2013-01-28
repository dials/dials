
import unittest

class TestRotationAngles(unittest.TestCase):

    def setUp(self):
        from scitbx import matrix
        from dials.array_family import flex
        from dials.io import xdsio
        from math import ceil, pi

        # The XDS files to read from    
        integrate_filename = './test/data/sim_mx/INTEGRATE.HKL'
        gxparm_filename = './test/data/sim_mx/GXPARM.XDS'
        
        # Read the XDS files
        self.integrate_handle = xdsio.IntegrateFile()
        self.integrate_handle.read_file(integrate_filename)
        self.gxparm_handle = xdsio.GxParmFile()
        self.gxparm_handle.read_file(gxparm_filename)

        # Get the parameters we need from the GXPARM file
        self.beam = self.gxparm_handle.get_beam()
        self.gonio = self.gxparm_handle.get_goniometer()
        self.detector = self.gxparm_handle.get_detector()
        self.ub_matrix = self.gxparm_handle.get_ub_matrix()
        self.ub_matrix = tuple(matrix.sqr(self.ub_matrix).inverse())
        self.unit_cell = self.gxparm_handle.get_unit_cell()
        self.space_group_type = self.gxparm_handle.get_space_group_type()

        # Get the minimum resolution in the integrate file
        d = [self.unit_cell.d(h) for h in self.integrate_handle.hkl]
        self.d_min = min(d)
        
        # Get the number of frames from the max z value 
        xcal, ycal, zcal = zip(*self.integrate_handle.xyzcal)
        self.gonio.num_frames = int(ceil(max(zcal)))
        
    def test_cctbx_angles(self):

        # FIXME this test is broken as there are different coordinate
        # frames in use - see coordinate_frame_converter in rstbx.cftbx
        
        from rstbx.diffraction import rotation_angles   
        from scitbx import matrix 
        
        ra = rotation_angles(self.d_min, 
                             self.ub_matrix, 
                             self.beam.wavelength, 
                             self.gonio.rotation_axis)
            
        ub = matrix.sqr(self.ub_matrix)
        s0 = matrix.col(self.beam.direction)
        m2 = matrix.col(self.gonio.rotation_axis)
        h = matrix.col(self.integrate_handle.hkl[0])
                         
        if ra(h):
            angles = ra.get_intersection_angles()
            
            for phi in angles:
                r = m2.axis_and_angle_as_r3_rotation_matrix(angle=phi)
                pstar = r * ub * h
                s1 = s0 + pstar
                self.assertAlmostEqual(s1.length(), s0.length())

    def test_xds_angles(self):
        from dials.spot_prediction import XdsRotationAngles   
        from scitbx import matrix 
        
        ra = XdsRotationAngles(self.beam.direction,
                               self.gonio.rotation_axis)
            
        ub = matrix.sqr(self.ub_matrix)
        s0 = matrix.col(self.beam.direction)
        m2 = matrix.col(self.gonio.rotation_axis)
        h = matrix.col(self.integrate_handle.hkl[0])
                  
        angles = ra.calculate(h, ub)
        
        for phi in angles:
            r = m2.axis_and_angle_as_r3_rotation_matrix(angle=phi)
            pstar = r * ub * h
            s1 = s0 + pstar
            self.assertAlmostEqual(s1.length(), s0.length())

if __name__ == '__main__':
    unittest.main()
