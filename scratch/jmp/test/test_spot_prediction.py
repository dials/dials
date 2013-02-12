import unittest

class TestSpotPrediction(unittest.TestCase):
    """Test spot prediction."""

    def setUp(self):
        pass

    def test_from_hkl(self):
        from scitbx import matrix
        from dials.array_family import flex
        from dials.io import xdsio
        from dials.spot_prediction import RotationAngles
        from dials.spot_prediction import is_angle_in_range
        from dials.geometry.transform import FromHklToRsv
        from dials.geometry.transform import FromBeamVectorToImageVolume
        from math import ceil, pi

        # The XDS files to read from
        integrate_filename = './test/data/sim_mx/INTEGRATE.HKL'
        gxparm_filename = './test/data/sim_mx/GXPARM.XDS'

        # Read the XDS files
        integrate_handle = xdsio.IntegrateFile()
        integrate_handle.read_file(integrate_filename)
        gxparm_handle = xdsio.GxParmFile()
        gxparm_handle.read_file(gxparm_filename)

        # Get the parameters we need from the GXPARM file
        d_min = 1.0
        beam = gxparm_handle.get_beam()
        gonio = gxparm_handle.get_goniometer()
        detector = gxparm_handle.get_detector()
        ub_matrix = gxparm_handle.get_ub_matrix()
        ub_matrix = tuple(matrix.sqr(ub_matrix).inverse())

        # Get the number of frames from the max z value
        xcal, ycal, zcal = zip(*integrate_handle.xyzcal)
        gonio.num_frames = int(ceil(max(zcal)))

        # Get the list of miller indices from the HKL file
        miller_indices = flex.miller_index(integrate_handle.hkl)

        # Calculate the rotation angles
        valid_angle = flex.bool()
        rotation_angle_calculator = RotationAngles(d_min, ub_matrix, beam.wavelength,
                                                   gonio.rotation_axis)
        rotation_angles = rotation_angle_calculator.calculate(miller_indices,
                                                              valid_angle)

        # Ensure none is invalid
        for va in valid_angle:
            self.assertTrue(va)

        # Convert rotation angles to degrees
        r2d = 180.0 / pi
        for i in range(len(rotation_angles)):
            rotation_angles[i] = rotation_angles[i] * r2d % 360

        # Get the rotation range
        phi0 = gonio.starting_angle
        phi1 = gonio.get_angle_from_frame(gonio.num_frames + 1)

        # Check which angles are in range
        in_range = []
        for phi in rotation_angles:
            in_range.append(phi0 <= phi <= phi1)

        # Ensure we only have those rotation angles that are valid
        new_rotation_angles = flex.double(len(miller_indices))
        for i in range(len(miller_indices)):
            if i > 0 and miller_indices[i] == miller_indices[i-1]:
                self.assertTrue(in_range[i+len(miller_indices)])
                new_rotation_angles[i] = rotation_angles[i+len(miller_indices)]
            else:
                self.assertTrue(in_range[0] or in_range[i+len(miller_indices)])
                if in_range[i]:
                    new_rotation_angles[i] = rotation_angles[i]
                elif in_range[i+len(miller_indices)]:
                    new_rotation_angles[i] = rotation_angles[i+len(miller_indices)]

        rotation_angles = new_rotation_angles

        # Get the reciprocal space vectors and vbeam vectors
        from_hkl_to_rsv = FromHklToRsv(ub_matrix, gonio.rotation_axis)
        rsv = from_hkl_to_rsv.apply(miller_indices, rotation_angles)
        beam_vectors = rsv + beam.direction

        # Calculate the image volume coordiantes
        status = flex.bool()
        from_beam_vector_to_image_volume = FromBeamVectorToImageVolume(detector, gonio)
        coords = from_beam_vector_to_image_volume.apply(beam_vectors, rotation_angles, status)

        # Get the XDS coords for comparison
        xds_coords = integrate_handle.xyzcal

        # Check that all succeded
        for i, s in enumerate(status):
            self.assertTrue(s)

        # Check the XY difference
        for xyz1, xyz2 in zip(coords, xds_coords):
            xy1 = matrix.col((xyz1[0], xyz1[1]))
            xy2 = matrix.col((xyz2[0], xyz2[1]))
            diff = (xy1 - xy2).length()
            self.assertLess(diff, 0.1)

        # Check the phi difference
        for xyz1, xyz2 in zip(coords, xds_coords):
            diff = abs(gonio.get_angle_from_frame(xyz1[2]+1) -
                       gonio.get_angle_from_frame(xyz2[2]+1))
            self.assertLess(diff, 0.1)

if __name__ == '__main__':
    unittest.main()
