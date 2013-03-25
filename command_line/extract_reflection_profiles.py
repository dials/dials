#!/usr/bin/env python
#
# extract_reflection_profiles.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

def print_ub_matrix(UB):
    '''Print the UB matrix to stdout.'''
    print "UB Matrix:"
    print "    (({0}, {1}, {2}),".format(UB[0], UB[1], UB[2])
    print "     ({0}, {1}, {2}),".format(UB[3], UB[4], UB[5])
    print "     ({0}, {1}, {2}))".format(UB[6], UB[7], UB[8])

def print_call_info(callback, info, result_type):
    '''Call a function and print basic info about the time it took to run the
    function and the number of output elements.'''
    from time import time
    print ""
    print info + ":"
    start_time = time()
    result = callback()
    time_taken = time() - start_time
    if result_type:
        print "{0} {1} in {2} s".format(len(result), result_type, time_taken)
    return result

def print_reflection_stats(reflections, adjacency_list):
    '''Print some reflection statistics.'''
    import numpy

    # Get the stats
    num_reflections = len(reflections)
    spot_x = [r.image_coord_px[0] for r in reflections]
    spot_y = [r.image_coord_px[1] for r in reflections]
    spot_z = [r.frame_number for r in reflections]
    bounding_boxes = [r.bounding_box for r in reflections]

    # Calculate the min, max, mean pixels in bounding box
    bbox_count = [(s[1]-s[0])*(s[3]-s[2])*(s[5]-s[4]) for s in bounding_boxes]
    min_bbox_size = numpy.min(bbox_count)
    max_bbox_size = numpy.max(bbox_count)
    med_bbox_size = int(numpy.median(bbox_count))

    # Calculate the mib, max, mean fast range of bbox
    bbox_fast_range = [s[1] - s[0] for s in bounding_boxes]
    min_bbox_fast_range = numpy.min(bbox_fast_range)
    max_bbox_fast_range = numpy.max(bbox_fast_range)
    med_bbox_fast_range = int(numpy.median(bbox_fast_range))

    # Calculate the mib, max, mean slow range of bbox
    bbox_slow_range = [s[3] - s[2] for s in bounding_boxes]
    min_bbox_slow_range = numpy.min(bbox_slow_range)
    max_bbox_slow_range = numpy.max(bbox_slow_range)
    med_bbox_slow_range = int(numpy.median(bbox_slow_range))

    # Calculate the mib, max, mean frame range of bbox
    bbox_frame_range = [s[5] - s[4] for s in bounding_boxes]
    min_bbox_frame_range = numpy.min(bbox_frame_range)
    max_bbox_frame_range = numpy.max(bbox_frame_range)
    med_bbox_frame_range = int(numpy.median(bbox_frame_range))

    # Get min/max/med bbox ranges
    min_bbox_range = (min_bbox_fast_range,
                      min_bbox_slow_range,
                      min_bbox_frame_range)
    max_bbox_range = (max_bbox_fast_range,
                      max_bbox_slow_range,
                      max_bbox_frame_range)
    med_bbox_range = (med_bbox_fast_range,
                      med_bbox_slow_range,
                      med_bbox_frame_range)

    bbox_count = (min_bbox_size, max_bbox_size, med_bbox_size)

    # Print the stats
    print ""
    print "Num reflections:", num_reflections
    print "Max spot x/y/z:", max(spot_x), max(spot_y), max(spot_z)
    print "Min spot x/y/z:", min(spot_x), min(spot_y), min(spot_z)
    print "Min/Max/Median bbox element count: ", bbox_count
    print "Max bbox range: ", max_bbox_range
    print "Min bbox range: ", min_bbox_range
    print "Med bbox range: ", med_bbox_range
    print "Num overlaps: ", adjacency_list.num_edges()

def run(xparm_path, integrate_path, image_frames, interactive):
    """Read the required data from the file, predict the spots and display."""

    from dials.algorithms.spot_prediction import IndexGenerator
    from dials.algorithms.spot_prediction import RayPredictor
    from dials.algorithms.spot_prediction import ray_intersection
    from dials.algorithms.spot_prediction import reflection_frames
    from dials.algorithms.integration import BBoxCalculator
    from dials.algorithms.integration import find_overlapping_reflections
    from dials.algorithms.integration import extract_reflection_profiles
    from iotbx.xds import xparm
    from dials.util import io
    from math import pi
    import dxtbx
    from dxtbx.sweep import SweepFactory
    from iotbx.xds import integrate_hkl

    # Set the number of frames
    num_frames = len(image_frames)

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

    # Read the sweep
    print "Reading: \"{0}\"".format("sweep")
    sweep = SweepFactory.sweep(image_frames)

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

    # Print the model data
    print ""
    print beam
    print detector
    print gonio
    print scan
    print "Resolution: ", d_min
    print "Sigma mosaicity: ", sigma_mosaicity
    print "Sigma Divergence: ", sigma_divergence
    print ""
    print_ub_matrix(UB)

    # Create the index generator
    generate_indices = IndexGenerator(unit_cell, space_group, d_min)

    # Create the spot predictor
    predict_rays = RayPredictor(beam.get_s0(), gonio.get_rotation_axis(), UB,
                                scan.get_oscillation_range(deg=False))

    # Generate Indices
    miller_indices = print_call_info(generate_indices.to_array,
        "Generating miller indices", "miller indices")

    # Predict reflections
    reflections = print_call_info(lambda: predict_rays(miller_indices),
        "Predicting rays", "reflections")

    # Get detector coordinates (mm)
    reflections = print_call_info(
        lambda: ray_intersection(detector, reflections),
        "Calculating detector coordinates", "coordinates")

    # Calculate the frame numbers of all the reflections
    reflections = print_call_info(
        lambda: reflection_frames(scan, reflections),
        "Calculating frame numbers", "frames")

    # Create the bbox calculator
    calculate_bbox = BBoxCalculator(beam, detector, gonio, scan,
        delta_divergence, delta_mosaicity)

    # Calculate the frame numbers of all the reflections
    reflections = print_call_info(
        lambda: (calculate_bbox(reflections), reflections)[1],
        "Calculating bboxes", "bboxes")

    # Find all the overlapping reflections
    adjacency_list = print_call_info(
        lambda: find_overlapping_reflections(reflections),
        "Calculating overlapping reflections", "edges")

    # Copy the reflection profiles from the sweep the reflection objects
    reflections = print_call_info(
        lambda: extract_reflection_profiles(sweep, reflections),
        "Copying reflection profiles from sweep", "reflections")

    # Print some reflection statistics
    print_reflection_stats(reflections, adjacency_list)

    # Enter an interactive python session
    if interactive:
        from dials.util.command_line import interactive_console
        interactive_console(namespace=locals())

if __name__ == '__main__':

    from optparse import OptionParser

    # Specify the command line options
    usage = "usage: %prog [options] /path/to/GXPARM.XDS /path/to/image.cbf"
    parser = OptionParser(usage)

    parser.add_option('-i', '--interactive',
                      dest='interactive', action="store_true", default=False,
                      help='Enter an interactive python session')

    # Parse the arguments
    (options, args) = parser.parse_args()

    # Print help if no arguments specified, otherwise call spot prediction
    if len(args) < 3:
        print parser.print_help()
    else:
        run(args[0], args[1], args[2:], options.interactive)
