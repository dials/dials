#!/usr/bin/env cctbx.python

def display_image_with_predicted_spots(image, xcoords, ycoords):
    """Display the image with coordinates overlayed."""
    from matplotlib import pylab, cm
    plt = pylab.imshow(image, vmin=0, vmax=1000, cmap=cm.Greys_r)
    pylab.scatter(xcoords, ycoords, marker='x')
    plt.axes.get_xaxis().set_ticks([])
    plt.axes.get_yaxis().set_ticks([])
    pylab.show()

def visualize_predicted_spots(image_volume, display_frame, spot_coords):
    """Get just those spots on the selected image and display."""
    spot_xy = [(x, y) for x, y, z in spot_coords if display_frame <= z < display_frame+1]
    xcoords, ycoords = zip(*spot_xy)
    display_image_with_predicted_spots(image_volume[display_frame,:,:], 
                                       xcoords, ycoords)

def predict_spots(input_filename, cbf_search_path, d_min, display_frame):
    """Read the required data from the file, predict the spots and display."""
    from dials.spot_prediction import SpotPredictor
    from dials.io import xdsio, pycbf_extra
    from cctbx import uctbx, sgtbx
    from time import time

    # Create the GXPARM file and read the contents
    print "Reading: \"{0}\"".format(input_filename)
    gxparm_handle = xdsio.GxParmFile()
    gxparm_handle.read_file(input_filename)
    beam      = gxparm_handle.get_beam() 
    gonio     = gxparm_handle.get_goniometer()
    detector  = gxparm_handle.get_detector()
    ub_matrix = gxparm_handle.get_ub_matrix()
    symmetry  = gxparm_handle.space_group

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
    spot_predictor = SpotPredictor(beam, gonio, detector, ub_matrix, d_min,
                                   unit_cell, space_group_type)
    
    # Predict the spot image volume coordinates 
    print "Predicting spots"
    start_time = time()
    spot_predictor.predict_spots()
    image_volume_coords = spot_predictor.image_volume_coordinates
    finish_time = time()
    print "Time taken: {0} s".format(finish_time - start_time)

    # If display frame selected then visualize
    if 0 <= display_frame < gonio.num_frames:
        print "Displaying predicted spots for frame \"{0}\"".format(display_frame)
        visualize_predicted_spots(image_volume, display_frame, 
                                  image_volume_coords)

if __name__ == '__main__':

    from optparse import OptionParser

    # Specify the command line options
    usage = "usage: %prog [options] /path/to/GXPARM.XDS"
    parser = OptionParser(usage)
    parser.add_option('-r', '--dmin',
                      dest='dmin',
                      type="float",
                      default=0.7,
                      help='Specify the resolution')
    parser.add_option('-c', '--cbf-search-path',
                      dest='cbf_search_path', 
                      default=None,
                      help='Specify search path for CBF files')
    parser.add_option('-d', '--display-frame',
                      dest='display_frame',
                      type="int",
                      default=-1,
                      help='Select a frame to display with predicted spots')
    (options, args) = parser.parse_args()

    # Print help if no arguments specified, otherwise call spot prediction
    if len(args) == 0:
        print parser.print_help()
    else:
        predict_spots(args[0], 
                      options.cbf_search_path, 
                      options.dmin, 
                      options.display_frame)
        
