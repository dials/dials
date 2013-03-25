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

def run(xparm_path, integrate_path, image_frames, interactive, output_file):
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
    from reflection_stats import ReflectionStats
    print ReflectionStats(reflections, adjacency_list)

    # Enter an interactive python session
    if interactive:
        from dials.util.command_line import interactive_console
        interactive_console(namespace=locals())

    # Dump the reflections to file
    if output_file:
        import pickle
        print "\nPickling the reflection list."
        pickle.dump(reflections, open(output_file, 'wb'))


if __name__ == '__main__':

    from optparse import OptionParser

    # Specify the command line options
    usage = "usage: %prog [options] /path/to/GXPARM.XDS /path/to/image.cbf"
    parser = OptionParser(usage)

    parser.add_option('-i', '--interactive',
                      dest='interactive', action="store_true", default=False,
                      help='Enter an interactive python session')
    parser.add_option('-o', '--output-file',
                      dest='output_file', type="string", default="",
                      help='Enter a destination filename for reflections')

    # Parse the arguments
    (options, args) = parser.parse_args()

    # Print help if no arguments specified, otherwise call spot prediction
    if len(args) < 3:
        print parser.print_help()
    else:
        run(args[0], args[1], args[2:], options.interactive,
            options.output_file)
