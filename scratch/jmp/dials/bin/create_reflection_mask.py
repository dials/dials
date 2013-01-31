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
    fig = pylab.figure()
    plt = pylab.imshow(mask[display_frame,:,:], cmap=cm.Greys_r, interpolation="nearest", origin='lower')
    #pylab.scatter(xcoords, ycoords, marker='x')
    plt.axes.get_xaxis().set_ticks([])
    plt.axes.get_yaxis().set_ticks([])
    pylab.show()
    #pylab.savefig('reflection_mask_frame_{0}.tiff'.format(display_frame), bbox_inches=0)
    #pylab.clf()

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
    from dials.integration import ReflectionMaskCreator
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
    reflections = spot_predictor.predict()
    finish_time = time()
    print "Time taken: {0} s".format(finish_time - start_time)
    
    # Create a mask of bad pixels
    detector_mask = flex.int(flex.grid(image_volume.shape[1:]), 0)
    detector_image = image_volume[0,:,:]
    for j in range(detector_mask.all()[0]):
        for i in range(detector_mask.all()[1]):
            if (detector_image[j,i] < 0):
                detector_mask[j,i] = -2    
    
    # Create the reflection mask    
    print "Creating reflection mask Roi for {0} spots".format(len(reflections))
    start_time = time()
    reflection_mask_creator = ReflectionMaskCreator(
                            beam, detector, gonio, 
                            detector_mask,                            
                            image_volume.shape,
                            sigma_divergence, 
                            sigma_mosaicity,
                            n_sigma)
    reflections = reflection_mask_creator.create(reflections)
    mask = reflection_mask_creator.mask
    finish_time = time()
    print "Time taken: {0} s".format(finish_time - start_time)
    
    # Extract arrays from array of reflections
    region_of_interest = []
    image_volume_coords = []
    for r in reflections:
        region_of_interest.append(r.region_of_interest)
        image_volume_coords.append(r.image_coord)
    
    # Get ranges of reflection masks
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
 
#    # Display frames   
#    if display_frame:
#        for frame in display_frame:
#            print "Displaying reflection mask for frame \"{0}\"".format(frame)
#            visualize_frame_reflection_mask(
#                reflection_mask_creator.mask.as_numpy_array(), 
#                frame, 
#                image_volume_coords)

#    # Display spots
#    if display_spot:
#        for spot in display_spot:
#            print "Displaying reflection mask for spot \"{0}\"".format(spot)
#            visualize_spot_reflection_mask(
#                reflection_mask_creator.mask.as_numpy_array(), 
#                spot, 
#                image_volume_coords, 
#                region_of_interest, 10)

    # 3D volume rendering
#    nx = 1000
#    ny = 1000
#    nz = 0
#    image_volume = image_volume[nz:,nx:,ny:]
#    for i in range(len(region_of_interest)):
#        roi = region_of_interest[i]
#        roi = (roi[0] - nx, roi[1] - nx, 
#               roi[2] - ny, roi[3] - ny, 
#               roi[4] - nz, roi[5] - nz)
#        region_of_interest[i] = roi
#    new_roi = []
#    for roi in region_of_interest:
#        if not (roi[0] < 0 or roi[2] < 0 or roi[4] < 0):
#            new_roi.append(roi)
#    region_of_interest = new_roi   
             
##    from scipy.ndimage.interpolation import zoom
#    import numpy
##    miller_indices = flex.miller_index(len(reflections))
##    for i, r in enumerate(reflections):
##        miller_indices[i] = r.miller_index
##    
##    count = 0
##    for i, hkl in enumerate(miller_indices):
##        if hkl == (30, -14, 12):
##            count = i
##            break
##            
##    print count
##    index = 28000
#    index = 1000#sim
##    factor = 1
#    roi = region_of_interest[index]
#    xyz = image_volume_coords[index]
#    image_volume = image_volume[xyz[2]-4:xyz[2]+4+1, roi[2]:roi[3], roi[0]:roi[1]]

#    from matplotlib import pylab, cm
#    for z in range(9):
#        fig = pylab.figure()
#        plt = pylab.imshow(image_volume[z,:,:], interpolation='nearest', 
#                     origin='lower', cmap=cm.Greys_r, 
#                     vmin=0, vmax=numpy.max(image_volume))
#        plt.axes.get_xaxis().set_ticks([])
#        plt.axes.get_yaxis().set_ticks([])                     
#        pylab.savefig("spot_intensity_sim_frame_no_axes_{0}.tiff".format(z))
#        pylab.show()
#    image_volume = zoom(image_volume.astype(numpy.float32), factor)
#    region_of_interest = (0, image_volume.shape[2],
#                          0, image_volume.shape[1],
#                          0, image_volume.shape[0])
                 
    #from spot_visualization import SpotVisualization
    #vis = SpotVisualization()
    #vis.vmax = 2000#0.5 * numpy.max(image_volume)
    #vis.visualize_reflections(None, region_of_interest)
    #vis.visualize_reflections(image_volume, region_of_interest)

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
