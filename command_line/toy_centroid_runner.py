from __future__ import division

def toy_centroid_runner(xparm_file, integrate_hkl_file, image_file):
    '''From the geometry in the xparm file, the indices in integrate_hkl_file
    and the images corresponding to the sweep to be generated from the
    image file, calculate the shoeboxes and from there the centroids using the
    toy centroid code.'''

    from dials.algorithms.centroid.toy_centroid import toy_centroid
    from dials.algorithms.spot_prediction import IndexGenerator
    from dials.algorithms.spot_prediction import RayPredictor
    from dials.algorithms.spot_prediction import ray_intersection
    from dials.algorithms.spot_prediction import reflection_frames
    from dials.algorithms.integration import ShoeboxCalculator
    from iotbx.xds import xparm
    from dials.util import io
    from math import pi
    import dxtbx

    from predict_spots import print_ub_matrix, print_reflection_stats, \
        display_predicted_spots

    sweep = dxtbx.sweep(image_file)

    models = dxtbx.load(xparm_file)
    beam = models.get_beam()
    detector = models.get_detector()
    gonio = models.get_goniometer()
    scan = models.get_scan()
    #scan = sweep.get_scan()
    first_image = scan.get_image_range()[0]
    image_range = sweep.get_scan().get_image_range()
    scan.set_image_range(image_range)

    # Read other data (need to assume an XPARM file
    xparm_handle = xparm.reader()
    xparm_handle.read_file(xparm_file)
    UB = io.get_ub_matrix_from_xparm(xparm_handle)
    unit_cell = io.get_unit_cell_from_xparm(xparm_handle)
    space_group = io.get_space_group_type_from_xparm(xparm_handle)
    d_min = detector.get_max_resolution_at_corners(
        beam.get_direction(), beam.get_wavelength())

    # FIXME I need to get these from integrate_hkl

    n_sigma = 5.0
    delta_divergence = n_sigma * 0.016 * pi / 180.0
    delta_mosaicity = n_sigma * 0.008 * pi / 180.0

    # Create the shoebox calculator
    calculate_shoebox = ShoeboxCalculator(beam, detector, gonio, scan,
        delta_divergence, delta_mosaicity)

    from dials.model.data import Reflection, ReflectionList

    from cctbx.array_family import flex

    miller_indices = flex.miller_index()

    for record in open(integrate_hkl_file):
        if record.startswith('!'):
            continue

        values = map(float, record.split())

        i, sigi = values[3:5]

        if sigi < 0:
            continue

        if i / sigi < 40:
            continue

        h, k, l = map(int, map(round, values[:3]))

        miller_indices.append((h, k, l))

    # Create the spot predictor
    predict_rays = RayPredictor(beam.get_s0(), gonio.get_rotation_axis(), UB,
                                scan.get_oscillation_range())

    reflections = reflection_frames(scan, ray_intersection(
        detector, predict_rays(miller_indices)))

    # Calculate the frame numbers of all the reflections
    calculate_shoebox(reflections)

    # Print some reflection statistics
    print_reflection_stats(reflections)

    bounding_boxes = { }

    for r in reflections[:10]:
        miller = r.miller_index
        cmin, cmax, rmin, rmax, fmin, fmax = r.shoebox
        if not miller in bounding_boxes:
            bounding_boxes[miller] = []
        bounding_boxes[miller].append((fmin, fmax, rmin, rmax, cmin, cmax))

    # FIXME in here need to sort list by frame number

    tc = toy_centroid(bounding_boxes, sweep)
    centroids = tc.get_centroids()

    for hkl in centroids:
        for centroid in centroids[hkl]:
            print '%.1f %.1f %.1f %.1f %.1f %.1f' % centroid

if __name__ == '__main__':
    import sys

    xparm_file = sys.argv[1]
    integrate_hkl_file = sys.argv[2]
    image_file = sys.argv[3]

    toy_centroid_runner(xparm_file, integrate_hkl_file, image_file)
