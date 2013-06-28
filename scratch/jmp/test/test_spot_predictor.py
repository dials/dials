from __future__ import division

import unittest

class TestSpotPredictor(unittest.TestCase):

    first = True

    def init_test_stuff(self):
        from scitbx import matrix
        from dials.array_family import flex
        from dials.spot_prediction import SpotPredictor
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

        # Create the spot predictor
        self.spot_predictor = SpotPredictor(self.beam,
                                            self.detector,
                                            self.gonio,
                                            self.unit_cell,
                                            self.space_group_type,
                                            self.ub_matrix,
                                            self.d_min)

        # Predict the spot locations
        self.spot_predictor.predict()

    def setUp(self):
        if self.first:
            self.init_test_stuff()
            self.first = False

    def test_miller_index_set(self):
        """Ensure we have the whole set of miller indices"""
        gen_hkl = {}
        for hkl in self.spot_predictor.miller_indices:
            gen_hkl[hkl] = True
        for hkl in self.integrate_handle.hkl:
            self.assertTrue(gen_hkl[hkl])

    def test_rotation_angles(self):
        """Ensure the rotation angles agree with XDS"""
        from scitbx import matrix

        # Create a dict of lists of xy for each hkl
        gen_phi = {}
        for hkl, phi in zip(self.spot_predictor.miller_indices,
                            self.spot_predictor.rotation_angles):
            try:
                a = gen_phi[hkl]
                a.append(phi)
                gen_phi[hkl] = a
            except KeyError:
                gen_phi[hkl] = [phi]

        for hkl, xyz in zip(self.integrate_handle.hkl,
                            self.integrate_handle.xyzcal):

            xds_phi = self.gonio.get_angle_from_zero_based_frame(xyz[2])

            # Select the nearest xy to use if there are 2
            my_phi = gen_phi[hkl]
            if len(my_phi) == 2:
                my_phi0 = my_phi[0]
                my_phi1 = my_phi[1]
                diff0 = abs(xds_phi - my_phi0)
                diff1 = abs(xds_phi - my_phi1)
                if (diff0 < diff1):
                    my_phi = my_phi0
                else:
                    my_phi = my_phi1
            else:
                my_phi = my_phi[0]

            self.assertTrue(abs(xds_phi - my_phi) < 0.1)

    def test_beam_vectors(self):
        """Ensure |s1| == |s0|"""
        from scitbx import matrix
        s0_length = matrix.col(self.beam.direction).length()
        for s1 in self.spot_predictor.beam_vectors:
            s1_length = matrix.col(s1).length()
            self.assertAlmostEqual(s0_length, s1_length)

    def test_image_coordinates(self):
        """Ensure the image coordinates agree with XDS"""
        from scitbx import matrix

        # Create a dict of lists of xy for each hkl
        gen_xy = {}
        for hkl, xyz in zip(self.spot_predictor.miller_indices,
                            self.spot_predictor.image_coordinates):
            try:
                a = gen_xy[hkl]
                a.append((xyz[0], xyz[1]))
                gen_xy[hkl] = a
            except KeyError:
                gen_xy[hkl] = [(xyz[0], xyz[1])]

        for hkl, xyz in zip(self.integrate_handle.hkl,
                            self.integrate_handle.xyzcal):

            xds_xy = (xyz[0], xyz[1])

            # Select the nearest xy to use if there are 2
            my_xy = gen_xy[hkl]
            if len(my_xy) == 2:
                my_xy0 = my_xy[0]
                my_xy1 = my_xy[1]
                diff0 = (matrix.col(xds_xy) - matrix.col(my_xy0)).length()
                diff1 = (matrix.col(xds_xy) - matrix.col(my_xy1)).length()
                if (diff0 < diff1):
                    my_xy = my_xy0
                else:
                    my_xy = my_xy1
            else:
                my_xy = my_xy[0]

            if (abs(xds_xy[0] - my_xy[0]) > 0.1 or
                abs(xds_xy[1] - my_xy[1]) > 0.1):
                    print xds_xy, gen_xy[hkl]
            self.assertTrue(abs(xds_xy[0] - my_xy[0]) < 0.1)
            self.assertTrue(abs(xds_xy[1] - my_xy[1]) < 0.1)

#        # Create map of detector coords
#        print gonio
#        hkl_xyz = {}
#        for hkl, xyz in zip(integrate_handle.hkl, integrate_handle.xyzcal):
#            try:
#                a = hkl_xyz[hkl]
#                a.append(xyz)
#                hkl_xyz[hkl] = a
#            except KeyError:
#                hkl_xyz[hkl] = [xyz]
#
#        first = True
#        count = 0
#        for hkl, xyz in zip(spot_predictor.miller_indices, spot_predictor.image_coordinates):
#
#            try:
#                xds_xyz = hkl_xyz[hkl]
#            except KeyError:
#                continue
#
#
#            diff = (matrix.col(xyz) - matrix.col(xds_xyz[0])).length()
#            if (abs(diff) > 1):
#                if (len(xds_xyz) == 2):
#                    diff = (matrix.col(xyz) - matrix.col(xds_xyz[1])).length()
#                    if (abs(diff) > 1):
#                        count += 1
#                        if first:
#                            print xyz, xds_xyz[1], diff
#                            first = False
#                else:
#                    count += 1
#                    if first:
#                        print xyz, xds_xyz[0], diff
#                        first = False

#        print count
#        # Get the list of miller indices from the HKL file
#        #miller_indices = flex.miller_index(integrate_handle.hkl)


##    s0 = matrix.col(beam.direction)
##    UB = matrix.sqr(ub_matrix)
##    h = matrix.col(miller_indices[0])
##    pstar0 = UB * h
##
##    from rstbx import diffraction
##
##    from dials.spot_prediction import XdsRotationAngles
##    xds_ra = XdsRotationAngles(beam.direction, gonio.rotation_axis)
##    cctbx_ra = diffraction.rotation_angles(d_min, ub_matrix, beam.wavelength, gonio.rotation_axis)
##    angles = xds_ra.calculate(h, ub_matrix)#pstar0)
##    cctbx_ra(h)
##    angles2 = cctbx_ra.get_intersection_angles()
##    R = matrix.col(gonio.rotation_axis).axis_and_angle_as_r3_rotation_matrix(angle=angles[0])
###    R = matrix.col(gonio.rotation_axis).axis_and_angle_as_r3_rotation_matrix(angle=angles2[0])
##    pstar = R * pstar0
##    s1 = s0 + pstar
##    print s1, s1.length(), s0.length()
##
#if __name__ == '__main__':
#    test()

if __name__ == '__main__':
    unittest.main()
