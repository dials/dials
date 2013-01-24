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

def visualize_frame_reflection_mask(mask, display_frame, spot_coords):
    """Display the reflection mask for a given frame."""
    from matplotlib import pylab, cm
    spot_xy = [(x, y) for x, y, z in spot_coords if display_frame <= z < display_frame+1]
    xcoords, ycoords = zip(*spot_xy)
    fig = pylab.figure(figsize=(8,8))
    plt = pylab.imshow(mask[display_frame,:,:], cmap=cm.Greys_r, interpolation="nearest")
    #pylab.scatter(xcoords, ycoords, marker='x')
    plt.axes.get_xaxis().set_ticks([])
    plt.axes.get_yaxis().set_ticks([])
    pylab.show()

def visualize_spot_reflection_mask(mask, display_spot, image_volume_coords, 
                                   region_of_interest, padding):
    """Display the reflection mask for a given spot."""
    from matplotlib import pylab, cm
    roi = region_of_interest[display_spot]
    xyz = image_volume_coords[display_spot]
    print "Image volume coordinate: ", xyz
    print "Region of interest: ", roi
    roi2 = [max(roi[0] - padding, 0), min(roi[1] + padding, mask.shape[2]-1),
            max(roi[2] - padding, 0), min(roi[3] + padding, mask.shape[1]-1),
            max(roi[4] - padding, 0), min(roi[5] + padding, mask.shape[0]-1)]
    
    image = mask[xyz[2], roi2[2]:roi2[3]+1, roi2[0]:roi2[1]+1]
    pylab.imshow(image, cmap=cm.Greys_r, interpolation="nearest", 
        extent=[roi2[0], roi2[1], roi2[2], roi2[3]])
    pylab.plot([roi[0], roi[1], roi[1], roi[0], roi[0]],
               [roi[2], roi[2], roi[3], roi[3], roi[2]])
    pylab.show()

def filter_reflections_by_roi_volume(region_of_interest, percent):
    """Filter the reflections by roi volume.
    
    Calculate the volume of each reflection, filter out the 1.0-percent largest
    reflection volumes and return an array containing True/False if the volume
    if valid.
    
    Args:
        region_of_interest The array of rois
        percent The percent of volumes to keep (0.0 < percent <= 1.0)
    
    Returns:
        A boolean array, True/False is the roi valid
    
    """
    from scitbx.array_family import flex
    from heapq import nlargest

    # Check given percentage
    if percent <= 0 or percent > 1.0:
        raise ValueError
        
    # A Calculate the volume of each region of interest
    calculate_roi_volume = lambda roi: ((roi[1] - roi[0]) * 
                                        (roi[3] - roi[2]) * 
                                        (roi[5] - roi[4]))
    volume = map(calculate_roi_volume, region_of_interest)
    
    # Calculate the volume limit below which 99% of reflections are
    n_reflections = len(volume)
    volume_limit = nlargest(int((1.0 - percent) * n_reflections), volume)[-1]
    
    # Create an array which is true if reflection volume is below the limit
    result = flex.bool(n_reflections)
    for i, v in enumerate(volume):
        result[i] = v < volume_limit
    
    return result    

def create_reflection_mask(input_filename, cbf_search_path, d_min, 
                           sigma_divergence, sigma_mosaicity, n_sigma,
                           display_frame, display_spot):
    """Read the required data from the file, predict the spots and calculate
       gaussian model parameters."""
    from dials.spot_prediction import SpotPredictor
    from dials.integration import ReflectionMaskRoi
    from dials.integration import ReflectionMask
    from dials.io import xdsio, pycbf_extra
    from cctbx import uctbx, sgtbx
    from time import time
    from dials.array_family import remove_if_not

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
    finish_time = time()
    print "Time taken: {0} s".format(finish_time - start_time)

    # Get the data from the spot predictor
    miller_indices = spot_predictor.miller_indices
    rotation_angles = spot_predictor.rotation_angles
    beam_vectors = spot_predictor.beam_vectors
    image_volume_coords = spot_predictor.image_volume_coordinates

    # Create the reflection mask regions of interest
    print "Creating reflection mask Roi"
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
    print "Filtering reflections by ROI"
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
    print "Creating reflection mask"
    start_time = time()
    reflection_mask = ReflectionMask(image_volume.shape)
    reflection_mask.create(image_volume_coords, 
                           region_of_interest)
    finish_time = time()
    print "Time taken: {0} s".format(finish_time - start_time)
     
    if display_frame:
        for frame in display_frame:
            print "Displaying reflection mask for frame \"{0}\"".format(frame)
            visualize_frame_reflection_mask(reflection_mask.mask.as_numpy_array(), 
                frame, image_volume_coords)

    if display_spot:
        for spot in display_spot:
            print "Displaying reflection mask for spot \"{0}\"".format(spot)
            visualize_spot_reflection_mask(reflection_mask.mask.as_numpy_array(), 
                spot, image_volume_coords, region_of_interest, 20)

                              
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
    parser.add_option('-d', '--display-frame',
                      dest='display_frame',
                      type="string",
                      action="callback",
                      callback=display_frame_callback,
                      help='Select a frame to display the reflection mask')      
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
        create_reflection_mask(args[0], 
                               args[1], 
                               options.dmin,
                               options.sigma_d,
                               options.sigma_m,
                               options.n_sigma,
                               options.display_frame,
                               options.display_spot)
