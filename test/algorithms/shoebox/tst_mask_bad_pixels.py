from __future__ import division

class Test(object):

    def __init__(self):
        import libtbx.load_env
        import os

        try:
            dials_regression = libtbx.env.dist_path('dials_regression')
        except KeyError, e:
            print 'FAIL: dials_regression not configured'
            return

        self.sweep_filename = os.path.join(dials_regression,
            'centroid_test_data', 'sweep.json')

        self.crystal_filename = os.path.join(dials_regression,
            'centroid_test_data', 'crystal.json')


    def run(self):
        from dials.model.serialize import load
        from dials.algorithms.integration import allocate_reflection_profiles
        from dials.algorithms import shoebox
        from scitbx.array_family import flex

        # Load the sweep and crystal
        self.sweep = load.sweep(self.sweep_filename)
        self.crystal = load.crystal(self.crystal_filename)

        # Get the reflections and overlaps
        reflections = self.predict_reflections()

        # Allocate memory for reflection profiles
        reflections = allocate_reflection_profiles(reflections)

        # If the adjacency list is given, then create the reflection mask
        image_size = self.detector.get_image_size()
        detector_mask = self.sweep[0] >= 0
        shoebox_masker = shoebox.MaskBadPixels(detector_mask)
        shoebox_masker(reflections)

        # Run the tests
        self.tst_values(reflections, self.detector.get_image_size(),
            detector_mask)

    def tst_values(self, reflections, image_size, detector_mask):
        '''Ensure pixels are masked correctly.'''
        import numpy
        from dials.algorithms import shoebox
        from scitbx.array_family import flex

        # Check that all elements in non_overlapping masks are 1
        for r in reflections:
            mask = r.shoebox_mask
            x0, x1, y0, y1, z0, z1 = r.bounding_box
            expected = flex.int(flex.grid(z1 - z0, y1 - y0, x1 - x0), 0)
            for j in range(y1 - y0):
                for i in range(x1 - x0):
                    x = x0 + i
                    y = y0 + j
                    if x < 0 or y < 0 or x >= image_size[0] or y >= image_size[1]:
                        pass
                    elif detector_mask[y, x]:
                        for k in range(z1 - z0):
                            expected[k,j,i] = shoebox.MaskCode.Valid
            assert(expected == mask)

        # Passed that test
        print "OK"

    def predict_reflections(self):
        from dials.algorithms.spot_prediction import ray_intersection
        from dials.algorithms.spot_prediction import reflection_frames
        from dials.algorithms import shoebox
        from dials.algorithms.integration import extract_reflection_profiles
        from dials.algorithms.integration import filter_by_detector_mask
        from dials.algorithms.integration import filter
        from cctbx import sgtbx
        from dials.algorithms.spot_prediction import IndexGenerator
        from dials.algorithms.spot_prediction import RayPredictor
        from math import sqrt

        # Get models from the sweep
        self.beam = self.sweep.get_beam()
        self.detector = self.sweep.get_detector()
        self.gonio = self.sweep.get_goniometer()
        self.scan = self.sweep.get_scan()

        # Create the index generator
        self.generate_hkl = IndexGenerator(
            self.crystal.get_unit_cell(),
            sgtbx.space_group_type(self.crystal.get_space_group()),
            self.detector.get_max_resolution(self.beam.get_s0(),
                self.beam.get_wavelength()))

        # Create the spot predictor
        self.predict_rays = RayPredictor(
            self.beam.get_s0(),
            self.gonio.get_rotation_axis(),
            self.scan.get_oscillation_range(deg=False))

        # Create the bbox calculator
        self.compute_bbox = shoebox.BBoxCalculator(
            self.beam, self.detector, self.gonio, self.scan,
            5 * self.beam.get_sigma_divergence(deg=False),
            5 * self.crystal.get_mosaicity(deg=False))

        # Generate Indices
        miller_indices = self.generate_hkl.to_array()

        # Predict reflections
        reflections = self.predict_rays(miller_indices, self.crystal.get_A())

        # Get detector coordinates (mm)
        reflections = ray_intersection(self.detector, reflections)

        # Calculate the frame numbers of all the reflections
        reflections = reflection_frames(self.scan, reflections)

        # Calculate the bounding boxes of all the reflections
        self.compute_bbox(reflections)

        # Return the reflections and overlaps
        return reflections

if __name__ == '__main__':
    test = Test()
    test.run()
