from __future__ import division

def toy_centroid_runner(xparm_file, integrate_hkl_file, image_file, output_file):
    '''From the geometry in the xparm file, the indices in integrate_hkl_file
    and the images corresponding to the sweep to be generated from the
    image file, calculate the shoeboxes and from there the centroids using the
    toy centroid code.'''

    from dials.algorithms.centroid.toy_centroid_Lui import toy_centroid_lui
    from dials.algorithms.spot_prediction import IndexGenerator
    from dials.algorithms.spot_prediction import RayPredictor
    from dials.algorithms.spot_prediction import ray_intersection
    from dials.algorithms.spot_prediction import reflection_frames
    from dials.algorithms.integration import BBoxCalculator
    from dials.algorithms.integration import extract_reflection_profiles
    from iotbx.xds import xparm
    from dials.util import io
    from math import pi
    import dxtbx
    from dxtbx.sweep import SweepFactory
    from iotbx.xds import integrate_hkl

    from predict_spots import print_ub_matrix, print_reflection_stats, \
        display_predicted_spots
    sweep = SweepFactory.sweep(image_file)
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

    # Read the integrate file
    print "Reading: \"{0}\"".format(integrate_hkl_file)
    integrate_handle = integrate_hkl.reader()
    integrate_handle.read_file(integrate_hkl_file)
    sigma_divergence = integrate_handle.sigma_divergence
    sigma_mosaicity = integrate_handle.sigma_mosaicity

    # Set the divergence and mosaicity
    n_sigma = 5.0
    delta_divergence = n_sigma * sigma_divergence * pi / 180.0
    delta_mosaicity = n_sigma * sigma_mosaicity * pi / 180.0

    # Create the bounding box calculator
    calculate_bbox = BBoxCalculator(beam, detector, gonio, scan,
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

    print beam
    print gonio
    print UB
    print scan
    print detector

    # Create the spot predictor
    predict_rays = RayPredictor(beam.get_s0(), gonio.get_rotation_axis(), UB,
                                scan.get_oscillation_range(deg = False))

    reflections = reflection_frames(scan, ray_intersection(
        detector, predict_rays(miller_indices)))


    # Calculate the frame numbers of all the reflections
    calculate_bbox(reflections)

    # extract reflection profiles
    extract_reflection_profiles(sweep, reflections)

    # Print some reflection statistics
    print_reflection_stats(reflections, None)

#    bounding_boxes = { }

#    for r in reflections[:10]:
#        miller = r.miller_index
#        cmin, cmax, rmin, rmax, fmin, fmax = r.shoebox
#        if not miller in bounding_boxes:
#            bounding_boxes[miller] = []
#        bounding_boxes[miller].append((fmin, fmax, rmin, rmax, cmin, cmax))

    # FIXME in here need to sort list by frame number

    tc = toy_centroid_lui(reflections)
    reflections = tc.get_reflections()

    for ref in reflections:
        #centroid = ref.centroid_position + ref.centroid_variance
        #print '%.1f %.1f %.1f %.1f %.1f %.1f' % centroid
        print ref
    # Dump the reflections to file
    if output_file:
        import pickle
        print "\nPickling the reflection list."
        pickle.dump(reflections, open(output_file, 'wb'))

if __name__ == '__main__':

    from optparse import OptionParser

    # Specify the command line options
    usage = "usage: %prog [options] /path/to/GXPARM.XDS"
    usage += " /path/to/INTEGRATE.HKL /path/to/image.cbf"

    parser = OptionParser(usage)

    parser.add_option('-o', '--output-file',
                      dest = 'output_file', type = "string", default = "",
                      help = 'Enter a destination filename for reflections')

    # Parse the arguments
    (options, args) = parser.parse_args()

    # Print help if no arguments specified, otherwise call function
    if len(args) < 3:
        print parser.print_help()
    else:
        toy_centroid_runner(args[0], args[1], args[2:], options.output_file)


