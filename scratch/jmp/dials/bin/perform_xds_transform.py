#!/usr/bin/env python

def parse_list_string(string):
    """Parse a string in the following ways:
    string: 1, 2, 3        -> [1, 2, 3]
    string: 1 - 6          -> [1, 2, 3, 4, 5, 6]
    string: 1 - 6, 7, 8, 9 -> [1, 2, 3, 4, 5, 6, 7, 8, 9]
    """
    items = string.split(',')
    for i in range(len(items)):
        items[i] = items[i].split("-")
        if len(items[i]) == 1:
            items[i] = [int(items[i][0])]
        elif len(items[i]) == 2:
            items[i] = range(int(items[i][0]), int(items[i][1]) + 1)
        else:
            raise SyntaxError
    items = [item for sublist in items for item in sublist]
    return set(items)

def display_frame_callback(option, opt, value, parser):
    setattr(parser.values, option.dest, parse_list_string(value))

def display_spot_callback(option, opt, value, parser):
    setattr(parser.values, option.dest, parse_list_string(value))

def visualize_xds_transform(grid, display_spot):
    """Display the XDS transform for a given spot."""
    from matplotlib import pylab, cm
    import numpy
    spot_grid = grid[display_spot, :, :, :]
    for i in range(0,9):
        ax = pylab.subplot(3, 3, i+1)
        image = spot_grid[i, :, :]
        plt=pylab.imshow(image, vmin=0, vmax=numpy.max(spot_grid), 
            cmap=cm.Greys_r)#, interpolation='nearest')
        plt.axes.get_xaxis().set_ticks([])
        plt.axes.get_yaxis().set_ticks([])   
        ax.set_title("slice: {0}".format(i))         
    pylab.show()

def perform_xds_transform(input_filename, cbf_search_path, d_min, 
                          sigma_divergence, sigma_mosaicity, n_sigma,
                          display_spot):
    """Read the required data from the file, predict the spots and calculate
       gaussian model parameters."""
    from dials.spot_prediction import SpotPredictor
    from dials.integration import ReflectionMaskRoi
    from dials.integration import ReflectionMask
    from dials.integration import SubtractBackground
    from dials.integration import filter_reflections_by_roi_volume
    from dials.integration import XdsTransformGrid
    from dials.integration import XdsTransform
    from dials.io import xdsio, pycbf_extra
    from cctbx import uctbx, sgtbx
    from time import time
    from dials.array_family import remove_if_not
    from scitbx.array_family import flex

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
        image_volume = flex.int(image_volume)
        gonio.num_frames = image_volume.all()[0]

    # Create the spot predictor
    spot_predictor = SpotPredictor(beam, gonio, detector, ub_matrix, d_min,
                                   unit_cell, space_group_type)
    
    # Predict the spot image volume coordinates 
    print "Predicting spots"
    start_time = time()
    spot_predictor.predict_spots()
    finish_time = time()
    print "Time taken: {0} s".format(finish_time - start_time)

    # Get the data from the spot predictor
    miller_indices = spot_predictor.miller_indices
    rotation_angles = spot_predictor.rotation_angles
    beam_vectors = spot_predictor.beam_vectors
    image_volume_coords = spot_predictor.image_volume_coordinates

    # Create the reflection mask regions of interest
    n_reflections = len(image_volume_coords)
    print "Creating reflection mask Roi for {0} reflections".format(n_reflections)
    start_time = time()
    reflection_mask_roi = ReflectionMaskRoi(
                            beam, detector, gonio, 
                            n_sigma * sigma_divergence, 
                            n_sigma * sigma_mosaicity)
    region_of_interest = reflection_mask_roi.calculate(
                            beam_vectors, rotation_angles)
    finish_time = time()
    print "Time taken: {0} s".format(finish_time - start_time)

    # Filter the reflections by region of interest volume
    n_reflections = len(region_of_interest)
    print "Filtering {0} reflections by ROI".format(n_reflections)
    start_time = time()
    valid_roi = filter_reflections_by_roi_volume(region_of_interest, 0.99)
    miller_indices      = remove_if_not(miller_indices, valid_roi)
    rotation_angles     = remove_if_not(rotation_angles, valid_roi)
    beam_vectors        = remove_if_not(beam_vectors, valid_roi)
    image_volume_coords = remove_if_not(image_volume_coords, valid_roi)
    region_of_interest  = remove_if_not(region_of_interest, valid_roi)
    finish_time = time()
    print "Time taken: {0} s".format(finish_time - start_time)
    
    range_x = [roi[1] - roi[0] for roi in region_of_interest]
    range_y = [roi[3] - roi[2] for roi in region_of_interest]
    range_z = [roi[5] - roi[4] for roi in region_of_interest]
    volume = [rx * ry * rz for rx, ry, rz in zip(range_x, range_y, range_z)]
    print "Min/Max ROI X Range: ", min(range_x), max(range_x)
    print "Min/Max ROI Y Range: ", min(range_y), max(range_y)
    print "Min/Max ROI Z Range: ", min(range_z), max(range_z)
    print "Min/Max Roi Volume:  ", min(volume), max(volume)
        
    # Create the reflection mask itself
    n_reflections = len(region_of_interest)
    print "Creating reflection mask for {0} reflections".format(n_reflections)
    start_time = time()
    reflection_mask = ReflectionMask(image_volume.all())
    reflection_mask.create(image_volume_coords, region_of_interest)
    finish_time = time()
    print "Time taken: {0} s".format(finish_time - start_time)
     
    # Subtract the background for each reflection
    n_reflections = len(region_of_interest)
    print "Subtracting background for {0} reflections".format(n_reflections)
    start_time = time()
    subtract_background = SubtractBackground(image_volume, reflection_mask.mask, min_pixels=100)
    valid_background = flex.bool(len(region_of_interest))
    subtract_background.subtract(region_of_interest, valid_background)
    subtract_background.set_non_reflection_value(0)
    miller_indices      = remove_if_not(miller_indices, valid_background)
    rotation_angles     = remove_if_not(rotation_angles, valid_background)
    beam_vectors        = remove_if_not(beam_vectors, valid_background)
    image_volume_coords = remove_if_not(image_volume_coords, valid_background)
    region_of_interest  = remove_if_not(region_of_interest, valid_background)
    finish_time = time()
    print "Time taken: {0} s".format(finish_time - start_time)
     
    print "Initialising XDS Transform"
    start_time = time()
    n_reflections = len(region_of_interest)
    grid_origin = (4, 4, 4)
    xds_grid = XdsTransformGrid(n_reflections, grid_origin, sigma_divergence, 
                                sigma_mosaicity, n_sigma)
    xds_trans = XdsTransform(xds_grid, image_volume, reflection_mask.mask,
                             detector, beam, gonio)
    finish_time = time()
    print "Time taken: {0} s".format(finish_time - start_time)
    
    print "Performing XDS transform on {0} reflections".format(n_reflections)
    start_time = time()
    xds_trans.calculate(region_of_interest, beam_vectors, rotation_angles)
    finish_time = time()
    print "Time taken: {0} s".format(finish_time - start_time)
     
    if display_spot:
        for spot in display_spot:
            print "Displaying XDS Transform for spot \"{0}\"".format(spot)
            visualize_xds_transform(xds_grid.data.as_numpy_array(), spot)

                              
if __name__ == '__main__':

    from optparse import OptionParser

    # Specify the command line options
    usage = "usage: %prog [options] /path/to/GXPARM.XDS /cbf/search/path/*.cbf"
    parser = OptionParser(usage)
    parser.add_option('-r', '--dmin',
                      dest='dmin',
                      type="float",
                      default=0.7,
                      help='Specify the resolution')
    parser.add_option('--sigma-d',
                      dest='sigma_d',
                      type="float",
                      default=0.034,
                      help='Specify the standard deviation of the beam divergence')
    parser.add_option('--sigma-m',
                      dest='sigma_m',
                      type="float",
                      default=0.082,
                      help='Specify the standard deviation of the mosaicity')
    parser.add_option('--num-sigma',
                      dest='n_sigma',
                      type="float",
                      default=10,
                      help='Specify the number of standard deviations to use')
    parser.add_option('--display-spot',
                      dest='display_spot',
                      type="string",
                      action="callback",
                      callback=display_frame_callback,
                      help='Select a frame to display the reflection mask')                   
    (options, args) = parser.parse_args()

    # Print help if no arguments specified, otherwise call spot prediction
    if len(args) < 2:
        print parser.print_help()
    else:
        perform_xds_transform(args[0], 
                              args[1], 
                              options.dmin,
                              options.sigma_d,
                              options.sigma_m,
                              options.n_sigma,
                              options.display_spot)
