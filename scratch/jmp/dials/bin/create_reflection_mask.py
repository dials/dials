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
    spot_xy = [(x-0.5, y-0.5) for x, y, z in spot_coords if display_frame <= z < display_frame+1]
    xcoords, ycoords = zip(*spot_xy)
    fig = pylab.figure(figsize=(8,8))
    plt = pylab.imshow(mask[display_frame,:,:], cmap=cm.Greys_r, interpolation="nearest", origin='lower')
    pylab.scatter(xcoords, ycoords, marker='x')
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
        extent=[roi2[0], roi2[1], roi2[2], roi2[3]], origin='lower')
    #pylab.plot([roi[0], roi[1], roi[1], roi[0], roi[0]],
    #           [roi[2], roi[2], roi[3], roi[3], roi[2]])
    pylab.plot([roi[0]-0.5, roi[1]-0.5, roi[1]-0.5, roi[0]-0.5, roi[0]-0.5],
               [roi[2]-0.5, roi[2]-0.5, roi[3]-0.5, roi[3]-0.5, roi[2]-0.5])
    pylab.show()

def create_reflection_mask(input_filename, cbf_search_path, d_min, 
                           sigma_divergence, sigma_mosaicity, n_sigma,
                           display_frame, display_spot):
    """Read the required data from the file, predict the spots and calculate
       gaussian model parameters."""
    from dials.spot_prediction import SpotPredictor
    from dials.integration import ReflectionMaskRoi
    from dials.integration import ReflectionMask
    from dials.integration import filter_reflections_by_roi_volume
    from dials.io import xdsio, pycbf_extra
    from cctbx import uctbx, sgtbx
    from time import time
    from dials.array_family.flex import remove_if_not
    from dials.array_family import flex
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

    print beam
    print gonio
    print detector
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
    
    # Predict the spot image volume coordinates 
    print "Predicting spots"
    start_time = time()
    spot_predictor.predict()
    finish_time = time()
    print "Time taken: {0} s".format(finish_time - start_time)

    # Get the data from the spot predictor
    miller_indices = spot_predictor.miller_indices
    rotation_angles = spot_predictor.rotation_angles
    beam_vectors = spot_predictor.beam_vectors
    image_volume_coords = spot_predictor.image_coordinates

    print min(rotation_angles), max(rotation_angles)

    # Create the reflection mask regions of interest
    print "Creating reflection mask Roi for {0} spots".format(len(miller_indices))
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
    range_phi = [gonio.get_angle_from_zero_based_frame(roi[5]) -
                 gonio.get_angle_from_zero_based_frame(roi[4]) 
                    for roi in region_of_interest]
                 
    volume = [rx * ry * rz for rx, ry, rz in zip(range_x, range_y, range_z)]
    print "Min/Max ROI X Range:   ", min(range_x), max(range_x)
    print "Min/Max ROI Y Range:   ", min(range_y), max(range_y)
    print "Min/Max ROI Z Range:   ", min(range_z), max(range_z)
    print "Min/Max ROI Phi Range: ", min(range_phi), max(range_phi)
    print "Min/Max ROI Volume:    ", min(volume), max(volume)

    xyz = image_volume_coords[2000]
    xyz = (int(xyz[0]), int(xyz[1]), int(xyz[2]))
    roi = region_of_interest[2000]
    #region_of_interest[2000] = (xyz[0] - 5, xyz[0] + 5 + 1, xyz[1] - 5, xyz[1] + 5 + 1, xyz[2] - 1, xyz[2] + 1 + 1)
    region_of_interest[2000] = (roi[0], roi[1], roi[2], roi[3], xyz[2] - 1, xyz[2] + 1 + 1)      
      
    # Create the reflection mask itself
    print "Creating reflection mask for {0} reflections".format(len(region_of_interest))
    start_time = time()
    reflection_mask = ReflectionMask(image_volume.shape)
    valid_roi = reflection_mask.create(image_volume_coords, region_of_interest)
    #miller_indices      = remove_if_not(miller_indices, valid_roi)
    #rotation_angles     = remove_if_not(rotation_angles, valid_roi)
    #beam_vectors        = remove_if_not(beam_vectors, valid_roi)
    #image_volume_coords = remove_if_not(image_volume_coords, valid_roi)
    #region_of_interest  = remove_if_not(region_of_interest, valid_roi)
    finish_time = time()
    print "Time taken: {0} s".format(finish_time - start_time)
    print "Created reflection mask for {0} reflections".format(len(region_of_interest))
    
    #for i, roi in enumerate(region_of_interest):
    #    if volume[i] == max(volume):
    #        index = i
    
    index = 2000
    print index, region_of_interest[index], volume[index]
    
    hkl = miller_indices[index]
    roi = region_of_interest[index]
    xyz = image_volume_coords[index]
    #roi = (roi[0], roi[1]+1, roi[2], roi[3]+1, roi[4], roi[5]+1)
    #roi = (int(xyz[0]) - 8, int(xyz[0]) + 8 + 1, int(xyz[1]) - 7, int(xyz[1]) + 9 + 1, int(xyz[2]) - 20, int(xyz[2]) + 20 + 1)
    print hkl, xyz, roi 
    print xyz[0] - roi[0], roi[1] - xyz[0]
    print xyz[1] - roi[2], roi[3] - xyz[1]
    print xyz[2] - roi[4], roi[5] - xyz[2]
    print image_volume[xyz[2], roi[2]:roi[3], roi[0]:roi[1]]
    
    #index = 100
    from dials.integration import centroid3d, centroid2d
    #xyz_obs = centroid3d(flex.int(image_volume), roi)
    xyz_obs = centroid3d(flex.int(image_volume), reflection_mask.mask, roi, index)
    #xy_obs = centroid2d(flex.int(image_volume[xyz[2],:,:]), (roi[0], roi[1], roi[2], roi[3]))
    print xyz, xyz_obs

    from matplotlib import pylab, cm
#    pylab.imshow(image_volume[xyz[2], :, :], interpolation='nearest', origin='lower', cmap=cm.Greys_r, vmin=0, vmax=1000)
#    pylab.plot([roi[0]-0.5, roi[1]-0.5, roi[1]-0.5, roi[0]-0.5, roi[0]-0.5],
#               [roi[2]-0.5, roi[2]-0.5, roi[3]-0.5, roi[3]-0.5, roi[2]-0.5])
#    pylab.scatter(xyz[0]-0.5, xyz[1]-0.5, marker='x')
#    pylab.scatter(xyz_obs[0]-0.5, xyz_obs[1]-0.5, marker='x', color='red')
#    #pylab.scatter(xy_obs[0]-0.5, xy_obs[1]-0.5, marker='x', color='yellow')
#    pylab.show()
    import numpy
    total_image = numpy.zeros(dtype=numpy.int32, shape=(roi[3]-roi[2], roi[1]-roi[0]))
    for z in range(roi[4], roi[5]+1):
        image = image_volume[z,roi[2]:roi[3],roi[0]:roi[1]]
        total_image = total_image + image

    xy = centroid2d(flex.int(total_image), (0, roi[1]-roi[0], 0, roi[3]-roi[2]))
    print xy[0] + roi[0], xy[1] + roi[2]
    pylab.imshow(total_image, interpolation='nearest', origin='lower', cmap=cm.Greys_r, vmin=0, vmax=10000)
    pylab.show()
    
        
#    for z in range(roi[4], roi[5]+1):
#        mask = reflection_mask.mask.as_numpy_array()[z,roi[2]:roi[3],roi[0]:roi[1]]
#        ind = numpy.where(mask != index)
#        mask[ind] = 0
#        pylab.imshow(mask, interpolation='nearest', origin='lower', cmap=cm.Greys_r, vmin=0, vmax=1000)
#        image = image_volume[z,roi[2]:roi[3],roi[0]:roi[1]]
#        pylab.imshow(image, interpolation='nearest', origin='lower', cmap=cm.Greys_r, vmin=0, vmax=1000)
#        pylab.show()

#    image = image_volume[roi[4]:roi[5]+1, roi[2]:roi[3]+1, roi[0]:roi[1]+1]

#    print image    
    
    
#    if display_frame:
#        for frame in display_frame:
#            print "Displaying reflection mask for frame \"{0}\"".format(frame)
#            visualize_frame_reflection_mask(reflection_mask.mask.as_numpy_array(), 
#                frame, image_volume_coords)

#    if display_spot:
#        for spot in display_spot:
#            print "Displaying reflection mask for spot \"{0}\"".format(spot)
#            visualize_spot_reflection_mask(reflection_mask.mask.as_numpy_array(), 
#                spot, image_volume_coords, region_of_interest, 10)

                              
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
