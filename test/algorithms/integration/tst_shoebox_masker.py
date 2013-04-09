from __future__ import division

def tst_non_overlapping(reflections, non_overlapping):
    '''Ensure non-overlapping reflections have all their values 1.'''
    import numpy

    # Check that all elements in non_overlapping masks are 1
    for i in non_overlapping:
        mask = reflections[i].shoebox_mask.as_numpy_array()
        ind = numpy.where(mask != 1)[0]
        assert(len(ind) == 0)

    # Passed that test
    print "OK"

def tst_overlapping(reflections, overlapping, adjacency_list):
    '''Ensure masks for overlapping reflections are set properly.'''
    import numpy
    from scitbx import matrix

    # Loop through all overlaps
    for i in overlapping:
        r1 = reflections[i]
        bbox_1 = r1.bounding_box
        r1_coord = matrix.col(r1.image_coord_px + (r1.frame_number,))

        # Create a mask that we expect
        r1_size = (bbox_1[5] - bbox_1[4],
                   bbox_1[3] - bbox_1[2],
                   bbox_1[1] - bbox_1[0])
        expected_mask = numpy.ones(shape = r1_size, dtype=numpy.int32)

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

        # Check the masks are the same
        calculated_mask = r1.shoebox_mask.as_numpy_array()
        assert(numpy.all(calculated_mask == expected_mask))

    # Passed the test
    print "OK"

def tst_reflection_mask(reflections, adjacency_list):
    '''Ensure masked values are correct'''
    import numpy

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

    # Run tests of overlapping and non_overlapping reflections
    tst_non_overlapping(reflections, non_overlapping)
    tst_overlapping(reflections, overlapping, adjacency_list)

def run():
    """Read the required data from the file, predict the spots and display."""

    from dials.algorithms.spot_prediction import IndexGenerator
    from dials.algorithms.spot_prediction import RayPredictor
    from dials.algorithms.spot_prediction import ray_intersection
    from dials.algorithms.spot_prediction import reflection_frames
    from dials.algorithms.integration import BBoxCalculator
    from dials.algorithms.integration import find_overlapping_reflections
    from dials.algorithms.integration import extract_reflection_profiles
    from iotbx.xds import xparm
    from dials.util import ioutil
    from math import pi
    import dxtbx
    from dxtbx.sweep import SweepFactory
    from iotbx.xds import integrate_hkl
    import os

    import libtbx.load_env
    try:
        dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
        print 'FAIL: dials_regression not configured'
        return

    xparm_path = os.path.join(dials_regression, 'centroid_test_data',
                           'GXPARM.XDS')
    integrate_path = os.path.join(dials_regression, 'centroid_test_data',
                           'INTEGRATE.HKL')

    # Set the number of frames
    num_frames = 1000

    # Read the models from the input file
    print "Reading: \"{0}\"".format(xparm_path)
    models = dxtbx.load(xparm_path)
    beam = models.get_beam()
    detector = models.get_detector()
    gonio = models.get_goniometer()
    scan = models.get_scan()
    first_image = scan.get_image_range()[0]
    image_range = (first_image, first_image + num_frames)
    scan.set_image_range(image_range)

    # Read other data (need to assume an XPARM file
    xparm_handle = xparm.reader()
    xparm_handle.read_file(xparm_path)
    UB = io.get_ub_matrix_from_xparm(xparm_handle)
    unit_cell = io.get_unit_cell_from_xparm(xparm_handle)
    space_group = io.get_space_group_type_from_xparm(xparm_handle)
    d_min = detector.get_max_resolution_at_corners(
        beam.get_direction(), beam.get_wavelength())

    # If the integrate.hkl path has been set get the bbox parameters
    print "Reading: \"{0}\"".format(integrate_path)

    # Read the integrate file
    integrate_handle = integrate_hkl.reader()
    integrate_handle.read_file(integrate_path)

    # Get the sigma_divergance and mosaicity
    sigma_divergence = integrate_handle.sigma_divergence
    sigma_mosaicity = integrate_handle.sigma_mosaicity

    # Set the divergence and mosaicity
    n_sigma = 5.0
    delta_divergence = n_sigma * sigma_divergence * pi / 180.0
    delta_mosaicity = n_sigma * sigma_mosaicity * pi / 180.0

    # Create the index generator
    generate_indices = IndexGenerator(unit_cell, space_group, d_min)

    # Create the spot predictor
    predict_rays = RayPredictor(beam.get_s0(), gonio.get_rotation_axis(), UB,
                                scan.get_oscillation_range(deg=False))

    # Generate Indices
    miller_indices = generate_indices.to_array()

    # Predict reflections
    reflections = predict_rays(miller_indices)

    # Get detector coordinates (mm)
    reflections = ray_intersection(detector, reflections)

    # Calculate the frame numbers of all the reflections
    reflections = reflection_frames(scan, reflections)

    # Create the bbox calculator
    calculate_bbox = BBoxCalculator(beam, detector, gonio, scan,
        delta_divergence, delta_mosaicity)

    # Calculate the frame numbers of all the reflections
    calculate_bbox(reflections)

    # Find all the overlapping reflections
    adjacency_list = find_overlapping_reflections(reflections)

    from dials.algorithms.integration import allocate_reflection_profiles
    from dials.algorithms.integration import ShoeboxMasker

    # Allocate memory for reflection profiles
    reflections = allocate_reflection_profiles(reflections)

    # If the adjacency list is given, then create the reflection mask
    shoebox_masker = ShoeboxMasker()
    shoebox_masker(reflections, adjacency_list)

    # Run the test
    tst_reflection_mask(reflections, adjacency_list)

if __name__ == '__main__':
    run()
