#!/usr/bin/env python

def parse_list_string( string ):
    """Parse a string in the following ways:
    string: 1, 2, 3        -> [1, 2, 3]
    string: 1 - 6          -> [1, 2, 3, 4, 5, 6]
    string: 1 - 6, 7, 8, 9 -> [1, 2, 3, 4, 5, 6, 7, 8, 9]
    """
    items = string.split( ',' )
    for i in range( len( items ) ):
        items[i] = items[i].split( "-" )
        if len( items[i] ) == 1:
            items[i] = [int( items[i][0] )]
        elif len( items[i] ) == 2:
            items[i] = range( int( items[i][0] ), int( items[i][1] ) + 1 )
        else:
            raise SyntaxError
    items = [item for sublist in items for item in sublist]
    return set( items )

def display_frame_callback( option, opt, value, parser ):
    setattr( parser.values, option.dest, parse_list_string( value ) )

def display_image_with_predicted_spots( image, xcoords, ycoords ):
    """Display the image with coordinates overlayed."""
    from matplotlib import pylab, cm
    from matplotlib import transforms

    plt = pylab.imshow( image, vmin = 0, vmax = 1000, cmap = cm.Greys_r,
                       interpolation = 'nearest', origin = 'lower' )
    pylab.scatter( xcoords, ycoords, marker = 'x' )
    plt.axes.get_xaxis().set_ticks( [] )
    plt.axes.get_yaxis().set_ticks( [] )
    pylab.show()

def visualize_predicted_spots( image_volume, display_frame, spot_coords ):
    """Get just those spots on the selected image and display."""
    print '_________________________________________________________________'
    spot_xy = [( x - 0.5, y - 0.5 ) for x, y, z in spot_coords if display_frame <= z < display_frame + 1]
    xcoords, ycoords = zip( *spot_xy )
    display_image_with_predicted_spots( image_volume[display_frame, :, :], xcoords, ycoords )

    my_tst_code( image_volume[display_frame, :, :], xcoords, ycoords )

def my_tst_code( image2d, x_ls, y_ls ):
    import numpy
    import time
    # import ind_2d_integrate_tst01
    import ind_2d_integrate
    cntrd_xcoord = numpy.zeros( len( x_ls ) )
    cntrd_ycoord = numpy.zeros( len( x_ls ) )
    x_sigma = numpy.zeros( len( x_ls ) )
    y_sigma = numpy.zeros( len( x_ls ) )

    # print time.time()
    time1 = time.time()
    # print "time1 =", time1

    ind_2d_integrate_tst01.start( image2d, x_ls, y_ls , cntrd_xcoord, cntrd_ycoord, x_sigma, y_sigma )
    # ind_2d_integrate.start( image2d, x_ls, y_ls , cntrd_xcoord, cntrd_ycoord, x_sigma, y_sigma )

    time2 = time.time()
    # print "time2 =", time2
    timedif = float( time2 - time1 )
    print "timedif =", timedif


    for i in range( len( x_ls ) ):
        print 'x,y (centroid) =', cntrd_xcoord[i], ', ', cntrd_ycoord[i]
        print 'sigma ( x, y) =', x_sigma [i], ', ', y_sigma[i]

def predict_spots( input_filename, cbf_search_path, d_min, display_frame ):
    """Read the required data from the file, predict the spots and display."""
    from dials_jmp.spot_prediction import SpotPredictor
    from dials_jmp.io import xdsio, pycbf_extra
    from cctbx import uctbx, sgtbx
    from time import time
    from scitbx import matrix

    # Create the GXPARM file and read the contents
    print "Reading: \"{0}\"".format( input_filename )
    gxparm_handle = xdsio.GxParmFile()
    gxparm_handle.read_file( input_filename )
    beam = gxparm_handle.get_beam()
    gonio = gxparm_handle.get_goniometer()
    detector = gxparm_handle.get_detector()
    ub_matrix = gxparm_handle.get_ub_matrix()
    symmetry = gxparm_handle.space_group

    print beam
    print gonio
    print detector
    print "UB Matrix:"
    print "    (({0}, {1}, {2}),".format( ub_matrix[0], ub_matrix[1], ub_matrix[2] )
    print "     ({0}, {1}, {2}),".format( ub_matrix[3], ub_matrix[4], ub_matrix[5] )
    print "     ({0}, {1}, {2}))".format( ub_matrix[6], ub_matrix[7], ub_matrix[8] )
    print "Symmetry: ", symmetry
    print "D min: ", d_min

    # Create the unit cell and space group objects
    unit_cell = uctbx.unit_cell( orthogonalization_matrix = ub_matrix )
    space_group_type = sgtbx.space_group_type( sgtbx.space_group( 
                            sgtbx.space_group_symbols( symmetry ).hall() ) )

    # Load the image volume from the CBF files and set the number of frames
    if cbf_search_path:
        print "Searching \"{0}\" for CBF files".format( cbf_search_path )
        image_volume = pycbf_extra.search_for_image_volume( cbf_search_path )
        gonio.num_frames = image_volume.shape[0]

    # Create the spot predictor
    spot_predictor = SpotPredictor( beam, detector, gonio, unit_cell,
                                   space_group_type,
                                   matrix.sqr( ub_matrix ).inverse(), d_min )

    # Predict the spot image volume coordinates
    print "Predicting spots"
    start_time = time()
    reflections = spot_predictor.predict()
    finish_time = time()
    print "Time taken: {0} s".format( finish_time - start_time )

    image_volume_coords = []
    for r in reflections:
        image_volume_coords.append( r.image_coord )

    # If display frame selected then visualize
    for frame in display_frame:
        print "Displaying predicted spots for frame \"{0}\"".format( frame )
        visualize_predicted_spots( image_volume, frame,
                                  image_volume_coords )

if __name__ == '__main__':

    from optparse import OptionParser

    # Specify the command line options
    usage = "usage: %prog [options] /path/to/GXPARM.XDS"
    parser = OptionParser( usage )
    parser.add_option( '-r', '--dmin',
                      dest = 'dmin',
                      type = "float",
                      default = 0.7,
                      help = 'Specify the resolution' )
    parser.add_option( '-c', '--cbf-search-path',
                      dest = 'cbf_search_path',
                      default = None,
                      help = 'Specify search path for CBF files' )
    parser.add_option( '-d', '--display-frame',
                      dest = 'display_frame',
                      type = "string",
                      action = "callback",
                      callback = display_frame_callback,
                      help = 'Select a frame to display with predicted spots' )
    ( options, args ) = parser.parse_args()

    # Print help if no arguments specified, otherwise call spot prediction
    if len( args ) == 0:
        print parser.print_help()
    else:
        predict_spots( args[0],
                      options.cbf_search_path,
                      options.dmin,
                      options.display_frame )
