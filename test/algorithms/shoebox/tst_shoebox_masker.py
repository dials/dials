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
        reflections, adjacency_list = self.predict_reflections()

        # Allocate memory for reflection profiles
        reflections = allocate_reflection_profiles(reflections)

        # If the adjacency list is given, then create the reflection mask
        image_size = self.detector.get_image_size()
        detector_mask = flex.bool(flex.grid(image_size[1], image_size[0]), True)
        shoebox_masker = shoebox.Masker(detector_mask)
        shoebox_masker(reflections, adjacency_list)

        # Loop through all edges
        overlapping = []
        for e in adjacency_list.edges():
            v1, v2 = adjacency_list[e]
            overlapping.append(v1)
            overlapping.append(v2)

        # Ensure elements are unique
        overlapping = set(overlapping)

        # Ensure we have some overlaps
        assert(len(overlapping) > 0)

        # Get all non-overlapping reflections
        all_r = set(range(len(reflections)))
        non_overlapping = all_r.difference(overlapping)

        # Run the tests
        self.tst_non_overlapping(reflections, non_overlapping,
            self.detector.get_image_size())
        self.tst_overlapping(reflections, overlapping, adjacency_list,
            image_size)

    def tst_non_overlapping(self, reflections, non_overlapping, image_size):
        '''Ensure non-overlapping reflections have all their values 1.'''
        import numpy
        from dials.algorithms import shoebox

        # Check that all elements in non_overlapping masks are 1
        for i in non_overlapping:
            mask = reflections[i].shoebox_mask.as_numpy_array()
            bbox = reflections[i].bounding_box
            if bbox[0] < 0:
                x0 = 0 - bbox[0]
            else:
                x0 = 0
            if bbox[2] < 0:
                y0 = 0 - bbox[2]
            else:
                y0 = 0
            if bbox[1] > image_size[0]:
                x1 = image_size[0] - bbox[1]
            else:
                x1 = bbox[1] - bbox[0]
            if bbox[3] > image_size[1]:
                y1 = image_size[1] - bbox[2]
            else:
                y1 = bbox[3] - bbox[2]

            ind = numpy.where(mask[:,y0:y1,x0:x1] != shoebox.MaskCode.Valid)[0]
            assert(len(ind) == 0)

        # Passed that test
        print "OK"

    def tst_overlapping(self, reflections, overlapping,
        adjacency_list, image_size):
        '''Ensure masks for overlapping reflections are set properly.'''
        import numpy
        from scitbx import matrix
        from dials.algorithms import shoebox

        # Loop through all overlaps
        for i in overlapping:
            r1 = reflections[i]
            bbox_1 = r1.bounding_box
            r1_coord = matrix.col(r1.image_coord_px + (r1.frame_number,))

            # Create a mask that we expect
            r1_size = (bbox_1[5] - bbox_1[4],
                       bbox_1[3] - bbox_1[2],
                       bbox_1[1] - bbox_1[0])
            expected_mask = numpy.zeros(shape = r1_size, dtype=numpy.int32)
            expected_mask[:,:,:] = shoebox.MaskCode.Valid

            # Loop through all reflections which this reflection overlaps
            for j in adjacency_list.adjacent_vertices(i):
                r2 = reflections[j]
                bbox_2 = r2.bounding_box
                r2_coord = matrix.col(r2.image_coord_px + (r2.frame_number,))

                # Get bounding box of intersection
                bbox_3 = (max(bbox_1[0], bbox_2[0]), min(bbox_1[1], bbox_2[1]),
                          max(bbox_1[2], bbox_2[2]), min(bbox_1[3], bbox_2[3]),
                          max(bbox_1[4], bbox_2[4]), min(bbox_1[5], bbox_2[5]))

                # Check intersection is valid
                assert(bbox_3[0] < bbox_3[1])
                assert(bbox_3[2] < bbox_3[3])
                assert(bbox_3[4] < bbox_3[5])

                # Get the coordinates are all mask values
                mask_coord = []
                for k in range(bbox_3[4], bbox_3[5]):
                    for j in range(bbox_3[2], bbox_3[3]):
                        for i in range(bbox_3[0], bbox_3[1]):
                            mask_coord.append(matrix.col((i+0.5, j+0.5, k+0.5)))

                dist = lambda a, m: numpy.array([(a - b).length() for b in m])

                # Find the indices in the intersection area where r2 is closer to
                # the point than r1
                ind = numpy.where(dist(r1_coord, mask_coord) >
                                  dist(r2_coord, mask_coord))[0]

                # Set the mask values for r1 where r2 is closer to 0
                k0, k1 = bbox_3[4] - bbox_1[4], bbox_3[5] - bbox_1[4]
                j0, j1 = bbox_3[2] - bbox_1[2], bbox_3[3] - bbox_1[2]
                i0, i1 = bbox_3[0] - bbox_1[0], bbox_3[1] - bbox_1[0]
                intersect_mask = expected_mask[k0:k1, j0:j1, i0:i1]
                intersect_mask_1d = intersect_mask.reshape((-1))
                intersect_mask_1d[ind] = 0
                intersect_mask[:,:] = intersect_mask_1d.reshape(intersect_mask.shape)
                expected_mask[k0:k1, j0:j1, i0:i1] = intersect_mask

                if bbox_1[0] < 0:
                    expected_mask[:,:,0:0-bbox_1[0]] = 0
                if bbox_1[2] < 0:
                    expected_mask[:,0:0-bbox_1[2],:] = 0
                if bbox_1[1] > image_size[0]:
                    expected_mask[:,:,image_size[0] - bbox_1[1]:] = 0
                if bbox_1[3] > image_size[1]:
                    expected_mask[:,image_size[1] - bbox_1[2]:,:] = 0

            # Check the masks are the same
            calculated_mask = r1.shoebox_mask.as_numpy_array()
            assert(numpy.all(calculated_mask == expected_mask))

        # Passed the test
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

        # Find overlapping reflections
        overlaps = shoebox.find_overlapping(reflections)

        # Return the reflections and overlaps
        return reflections, overlaps

if __name__ == '__main__':
    test = Test()
    test.run()
