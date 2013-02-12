
import unittest

class TestRotationAngles(unittest.TestCase):

    def setUp(self):
        from scitbx import matrix
        from dials.array_family import flex
        from dials.io import xdsio
        from math import ceil, pi

        # The XDS files to read from
        self.integrate_filename = './test/data/sim_mx/INTEGRATE.HKL'
        self.gxparm_filename = './test/data/sim_mx/GXPARM.XDS'

        # Read the XDS files
        self.integrate_handle = xdsio.IntegrateFile()
        self.integrate_handle.read_file(self.integrate_filename)
        self.gxparm_handle = xdsio.GxParmFile()
        self.gxparm_handle.read_file(self.gxparm_filename)

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
        from rstbx.diffraction import rotation_angles
        from rstbx.cftbx.coordinate_frame_converter import coordinate_frame_converter
        from scitbx import matrix

        cfc = coordinate_frame_converter(self.gxparm_filename)
        u, b = cfc.get_u_b(convention = cfc.ROSSMANN)
        axis = cfc.get('rotation_axis', convention = cfc.ROSSMANN)
        ub = u * b
        wavelength = cfc.get('wavelength')

        ra = rotation_angles(self.d_min,
                             ub,
                             wavelength,
                             axis)

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
