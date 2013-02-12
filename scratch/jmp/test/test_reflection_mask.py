
import unittest

class TestReflectionMask(unittest.TestCase):

    def setUp(self):
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

    def test_centroids(self):

        from dials.spot_prediction import SpotPredictor
        from dials.integration import ReflectionMask, ReflectionMaskRoi

        # Create the spot predictor
        spot_predictor = SpotPredictor(self.beam,
                                       self.detector,
                                       self.gonio,
                                       self.unit_cell,
                                       self.space_group_type,
                                       self.ub_matrix,
                                       self.d_min)

        # Predict the spot locations of the set of miller indices
        spot_predictor.predict(flex.miller_index(self.integrate_handle.hkl)

        # Create the reflection mask regions of interest
#        reflection_mask_roi = ReflectionMaskRoi(
#                                self.beam, self.detector, self.gonio,
#                                n_sigma * sigma_divergence,
#                                n_sigma * sigma_mosaicity)
#        region_of_interest = reflection_mask_roi.calculate(
#                                self.beam_vectors, self.rotation_angles)
#
#        # Create the reflection mask itself
#        reflection_mask = ReflectionMask(image_volume.shape)
#        valid_roi = reflection_mask.create(image_volume_coords, region_of_interest)
#
if __name__ == '__main__':
    unittest.main()
