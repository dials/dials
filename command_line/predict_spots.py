#!/usr/bin/env python

def predict_spots(filenames, verbose):
    """Read the required data from the file, predict the spots and display."""
    from dials.algorithms.spot_prediction import SpotPredictor
    from dials_jmp.io import xdsio, pycbf_extra
    from cctbx import uctbx, sgtbx
    from time import time
    from scitbx import matrix

    # Create the GXPARM file and read the contents
    print "Reading: \"{0}\"".format(input_filename)
    gxparm_handle = xdsio.GxParmFile()
    gxparm_handle.read_file(input_filename)
    beam      = gxparm_handle.get_beam()
    gonio     = gxparm_handle.get_goniometer()
    detector  = gxparm_handle.get_detector()
    ub_matrix = gxparm_handle.get_ub_matrix()
    symmetry  = gxparm_handle.space_group

    # Print out model details
    print beam
    print detector
    print gonio
    print scan
    
    # Print out other details
    print "UB Matrix:"
    print "    (({0}, {1}, {2}),".format(ub_matrix[0], ub_matrix[1], ub_matrix[2])
    print "     ({0}, {1}, {2}),".format(ub_matrix[3], ub_matrix[4], ub_matrix[5])
    print "     ({0}, {1}, {2}))".format(ub_matrix[6], ub_matrix[7], ub_matrix[8])
    print "Symmetry: ", symmetry
    print "D min: ", d_min

    # Create the unit cell and space group objects
    unit_cell = uctbx.unit_cell(orthogonalization_matrix = ub_matrix)
    space_group_type = sgtbx.space_group_type(sgtbx.space_group(
                            sgtbx.space_group_symbols(symmetry).hall()))

    # Load the image volume from the CBF files and set the number of frames
    if cbf_search_path:
        print "Searching \"{0}\" for CBF files".format(cbf_search_path)
        image_volume = pycbf_extra.search_for_image_volume(cbf_search_path)
        gonio.num_frames = image_volume.shape[0]

    # Create the spot predictor
    spot_predictor = SpotPredictor(beam, detector, gonio, unit_cell,
                                   space_group_type,
                                   matrix.sqr(ub_matrix).inverse(), d_min)

    # Predict the reflections
    print "Predicting spots"
    start_time = time()
    reflections = spot_predictor()
    finish_time = time()
    print "Time taken: {0} s".format(finish_time - start_time)
    print "Num spots: {0}".format(len(reflections))
    
    # If verbose then print reflections
    if verbose:
        for r in reflections:
            print r


if __name__ == '__main__':

    from optparse import OptionParser

    # Specify the command line options
    usage = "usage: %prog [options] /path/to/GXPARM.XDS"
    parser = OptionParser(usage)

    # Add a verbose option (False by default)
    parser.add_option('-v', action="store_true", dest="verbose", default=False
                      help='Print out reflection details (default = False)')

    # Parse the arguments
    (options, args) = parser.parse_args()

    # Print help if no arguments specified, otherwise call spot prediction
    if len(args) == 0:
        print parser.print_help()
    else:
        predict_spots(args[0], options.verbose)
