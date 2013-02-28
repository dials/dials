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

def predict_spots(input_filename, numframes, d_min, verbose):
    """Read the required data from the file, predict the spots and display."""

    from dials.algorithms.spot_prediction import IndexGenerator
    from dials.algorithms.spot_prediction import SpotPredictor
    from iotbx.xds import xparm
    from dials.util import io

    # Read the models from the input file
    print "Reading: \"{0}\"".format(input_filename)
    models = io.read_models_from_file(input_filename)
    beam = models.get_beam()
    detector = models.get_detector()
    gonio = models.get_goniometer()
    scan = models.get_scan()
    scan.image_range = (scan.image_range[0], scan.image_range[0] + num_frames)

    # Read other data (need to assume an XPARM file
    xparm_handle = xparm.reader()
    xparm_handle.read_file(input_filename)
    ub_matrix = io.get_ub_matrix_from_xparm(xparm_handle)
    unit_cell = io.get_unit_cell_from_xparm(xparm_handle)
    space_group = io.get_space_group_type_from_xparm(xparm_handle)

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
    print ""
    print_ub_matrix(ub_matrix)

    # Create the index generator
    generate_indices = IndexGenerator(unit_cell, space_group, True, d_min)

    # Create the spot predictor
    predict_spots = SpotPredictor(beam, detector, gonio, scan, ub_matrix)

    # Generate Indices
    miller_indices = print_call_info(generate_indices.to_array,
        "Generating miller indices", "miller indices")

    # Predict reflections
    reflections = print_call_info(lambda: predict_spots(miller_indices),
        "Predicting spots", "reflections")

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
        predict_spots(args[0], 4000, 0.7, options.verbose)
