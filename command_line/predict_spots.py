#!/usr/bin/env python
#
# predict_spots.py
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
    print info
    start_time = time()
    result = callback()
    finish_time = time()
    time_taken = finish_time - start_time
    print "{0} {1} in {2} s".format(len(result), result_type, time_taken)
    return result

def get_intersection(detector, reflection):
    reflection.image_coord_mm = detector[0].get_ray_intersection(
        reflection.beam_vector)
    return reflection

def get_frame_numbers(scan, reflection):
    from dials.model.data import Reflection
    fn = scan.get_frames_with_angle(reflection.rotation_angle)
    reflection_list = []
    for f in fn:
        r = Reflection(reflection)
        r.frame_number = f
        reflection_list.append(r)
    return reflection_list

def predict_spots(input_filename, num_frames, verbose):
    """Read the required data from the file, predict the spots and display."""

    from dials.algorithms.spot_prediction import IndexGenerator
    from dials.algorithms.spot_prediction import RayPredictor
    from dials.algorithms.spot_prediction import RayIntersector
    from iotbx.xds import xparm
    from dials.util import io
    import dxtbx

    # Read the models from the input file
    print "Reading: \"{0}\"".format(input_filename)
    models = dxtbx.load(input_filename)
    beam = models.get_beam()
    detector = models.get_detector()
    gonio = models.get_goniometer()
    scan = models.get_scan()
    first_image = scan.get_image_range()[0]
    image_range = (first_image, first_image + num_frames)
    scan.set_image_range(image_range)

    # Read other data (need to assume an XPARM file
    xparm_handle = xparm.reader()
    xparm_handle.read_file(input_filename)
    UB = io.get_ub_matrix_from_xparm(xparm_handle)
    unit_cell = io.get_unit_cell_from_xparm(xparm_handle)
    space_group = io.get_space_group_type_from_xparm(xparm_handle)
    d_min = detector.get_max_resolution_at_corners(
        beam.get_direction(), beam.get_wavelength())
    # Load the image volume from the CBF files and set the number of frames
#    if cbf_path:
#        print "Loading CBF files"
#        image_volume = pycbf_extra.search_for_image_volume(cbf_search_path)
#        scan.image_range = (scan.image_range[0],
#            scan.image_range[0] + image_volume.shape[0] - 1)

    # Print the model data
    print ""
    print beam
    print detector
    print gonio
    print scan
    print "Resolution: ", d_min
    print ""
    print_ub_matrix(UB)

    # Create the index generator
    generate_indices = IndexGenerator(unit_cell, space_group, d_min)

    # Create the spot predictor
    predict_rays = RayPredictor(beam.get_s0(), gonio.get_rotation_axis(), UB,
                                scan.get_oscillation_range())

    # Generate Indices
    miller_indices = print_call_info(generate_indices.to_array,
        "Generating miller indices", "miller indices")

    # Predict reflections
    reflections = print_call_info(lambda: predict_rays(miller_indices),
        "Predicting rays", "reflections")

    #intersection = lambda x: get_intersection(detector, x)
    #all_intersections = lambda: map(intersection, reflections)
    intersection = RayIntersector(detector)

    # Get detector coordinates (mm)
    reflections = print_call_info(lambda: intersection(reflections),
        "Calculating detector coordinates", "coordinates")

    frame_calculator = lambda x: get_frame_numbers(scan, x)
    frame_calculator_all = lambda: map(frame_calculator, reflections)

    reflections = print_call_info(frame_calculator_all,
        "Calculating frame numbers", "frame")


    # Get the reflection frame numbers


#        // Get the list of frames at which the reflection will be observed
#        // and add the predicted observations to the list of reflections
#        flex_double frames = scan_.get_frames_with_angle(phi);
#        for (std::size_t j = 0; j < frames.size(); ++j) {
#          reflection_type r;
#          r.set_miller_index(h);
#          r.set_rotation_angle(phi);
#          r.set_beam_vector(s1);
#          r.set_image_coord_px(xy_px);
#          r.set_image_coord_mm(xy_mm);
#          r.set_frame_number(frames[j]);
#          reflections.push_back(r);
#        }
#      }

if __name__ == '__main__':

    from optparse import OptionParser

    # Specify the command line options
    usage = "usage: %prog [options] /path/to/GXPARM.XDS"
    parser = OptionParser(usage)

    # Add a verbose option (False by default)
    parser.add_option('-v', action="store_true", dest="verbose", default=False,
                      help='Print out reflection details (default = False)')

    # Parse the arguments
    (options, args) = parser.parse_args()

    # Print help if no arguments specified, otherwise call spot prediction
    if len(args) == 0:
        print parser.print_help()
    else:
        predict_spots(args[0], 4000, options.verbose)
