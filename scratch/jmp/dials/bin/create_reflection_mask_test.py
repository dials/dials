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

def create_reflection_mask(input_filename, integrate_filename, cbf_search_path, d_min, 
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

    xcorr_filename = '/home/upc86896/Projects/dials/dials-svn/dials-code/scratch/jmp/test/data/X-CORRECTIONS.cbf'
    ycorr_filename = '/home/upc86896/Projects/dials/dials-svn/dials-code/scratch/jmp/test/data/Y-CORRECTIONS.cbf'

    # Create the GXPARM file and read the contents
    print "Reading: \"{0}\"".format(input_filename)
    gxparm_handle = xdsio.GxParmFile()
    gxparm_handle.read_file(input_filename)
    beam      = gxparm_handle.get_beam() 
    gonio     = gxparm_handle.get_goniometer()
    detector  = gxparm_handle.get_detector()
    ub_matrix = gxparm_handle.get_ub_matrix()
    symmetry  = gxparm_handle.space_group
    integrate_handle = xdsio.IntegrateFile()
    integrate_handle.read_file(integrate_filename)
    xcorr_handle = xdsio.XYCorrection()
    xcorr_handle.read_file(xcorr_filename)
    ycorr_handle = xdsio.XYCorrection()
    ycorr_handle.read_file(ycorr_filename)

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

    miller_indices = integrate_handle.hkl
    new_miller_indices = []
    for i, hkl in enumerate(miller_indices):
        if i > 0:
            if hkl != last_hkl:
                new_miller_indices.append(hkl)
        else:
            new_miller_indices.append(hkl)
        last_hkl = hkl
    
    miller_indices = flex.miller_index(new_miller_indices)
    
    # Predict the spot image volume coordinates 
    print "Predicting spots"
    start_time = time()
    spot_predictor.predict(miller_indices)
    finish_time = time()
    print "Time taken: {0} s".format(finish_time - start_time)

    # Get the data from the spot predictor
    miller_indices = spot_predictor.miller_indices
    rotation_angles = spot_predictor.rotation_angles
    beam_vectors = spot_predictor.beam_vectors
    image_volume_coords = spot_predictor.image_coordinates

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
#    miller_indices      = remove_if_not(miller_indices, valid_roi)
#    rotation_angles     = remove_if_not(rotation_angles, valid_roi)
#    beam_vectors        = remove_if_not(beam_vectors, valid_roi)
#    image_volume_coords = remove_if_not(image_volume_coords, valid_roi)
#    region_of_interest  = remove_if_not(region_of_interest, valid_roi)
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

    # Create the reflection mask itself
    print "Creating reflection mask for {0} reflections".format(len(region_of_interest))
    start_time = time()
    reflection_mask = ReflectionMask(image_volume.shape)
    valid_roi = reflection_mask.create(image_volume_coords, region_of_interest)
    finish_time = time()
    print "Time taken: {0} s".format(finish_time - start_time)
    print "Created reflection mask for {0} reflections".format(len(region_of_interest))
    
    index = 120
    hkl_gen = miller_indices[index]
    roi_gen = region_of_interest[index]
    xyz_gen = image_volume_coords[index]

    import numpy
    xcorr = xcorr_handle.get_correction(detector.size[::-1])
    ycorr = ycorr_handle.get_correction(detector.size[::-1])

    from dials.integration import centroid3d, centroid2d, centroid_reflection

    #roi_gen = (roi_gen[0], roi_gen[1], roi_gen[2], roi_gen[3], 46, 54)

    #roi_gen = (roi_gen[0], roi_gen[1], roi_gen[2], roi_gen[3], int(xyz_gen[2])-2, int(xyz_gen[2])+2)
#    import numpy 
#    spot_image = image_volume[roi_gen[4]:roi_gen[5], roi_gen[2]:roi_gen[3], roi_gen[0]:roi_gen[1]]

#    mask = (spot_image >= 10).astype(numpy.int32)
#    
#    zyx = (26, 8, 7)
##    print zyx
#    mask = markregion(mask, zyx)
#        
#    numpy.set_printoptions(threshold=numpy.nan)
##    print mask
#   
#    ind = numpy.where(mask != 2)
#    spot_image[ind] = 0

#    print spot_image

#    image_volume[roi_gen[4]:roi_gen[5], roi_gen[2]:roi_gen[3], roi_gen[0]:roi_gen[1]] = spot_image

#    

#    try:
#        xyz_obs = centroid3d(flex.int(image_volume), reflection_mask.mask, roi_gen, index)
#    except RuntimeError:
#        print "Failed to calculated observed xyz"
#        xyz_obs = None

    threshold = 10
    xyz_obs = centroid_reflection(flex.int(image_volume), reflection_mask.mask, 
                                  roi_gen, index, threshold)

    
#    for z in range(roi_gen[4],roi_gen[5]):
#        mask = flex.int(reflection_mask.mask.as_numpy_array()[z,:,:])
#        xy = centroid2d(flex.int(image_volume[z,:,:]), mask, roi_gen[0:4], index)
#        xy = (xy[0] + xcorr[xy[1], xy[0]], xy[1] + ycorr[xy[1], xy[0]]) 
#        print z, xy
        
    #print xcorr[xyz_obs[1]-5:xyz_obs[1]+5, xyz_obs[0]-5:xyz_obs[0]+5]

    #xyz_obs_corr = centroid(flex.int(image_volume), reflection_mask.mask, roi_gen, index, xcorr, ycorr)


    
    xyz_obs_corr = (xyz_obs[0] + xcorr[xyz_obs[1], xyz_obs[0]], 
                    xyz_obs[1] + ycorr[xyz_obs[1], xyz_obs[0]], 
                    xyz_obs[2])
    #xy_obs = centroid2d(flex.int(image_volume[xyz[2],:,:]), (roi[0], roi[1], roi[2], roi[3]))
    xds_hkl_xyz_cal = {}
    for hkl, xyz in zip(integrate_handle.hkl, integrate_handle.xyzcal): 
        xds_hkl_xyz_cal[hkl] = xyz

    xds_hkl_xyz_obs = {}
    for hkl, xyz in zip(integrate_handle.hkl, integrate_handle.xyzobs): 
        xds_hkl_xyz_obs[hkl] = xyz
    
    print "HKL: ", hkl_gen
    print "XYZ_CAL: ({0:.1f}, {1:.1f}, {2:.1f})".format(*xyz_gen)
    print "ROI: ", roi_gen
    print "XYZ_OBS:  ({0:.1f}, {1:.1f}, {2:.1f})".format(*xyz_obs)
    print "XYZ_OBS_CORR:  ({0:.1f}, {1:.1f}, {2:.1f})".format(*xyz_obs_corr)
    print "XYZ_CAL_XDS: ", xds_hkl_xyz_cal[hkl_gen]
    print "XYZ_OBS_XDS: ", xds_hkl_xyz_obs[hkl_gen]

def centroid(image, mask, roi, value, xcorr, ycorr):
    xc = 0.0
    yc = 0.0
    zc = 0.0
    count = 0.0
    for k in range(roi[4], roi[5]):
        for j in range(roi[2], roi[3]):
            for i in range(roi[0], roi[1]):
                if (mask[k, j, i] == value):
                    xc += (i+0.5+xcorr[j,i]) * image[k, j, i]
                    yc += (j+0.5+ycorr[j,i]) * image[k, j, i]
                    zc += (k+0.5) * image[k, j, i]
                    count += image[k, j, i]
    
    return (xc / count, yc / count, zc / count)

def markregion(mask, seed):
    z, y, x = seed
    if (z >= 0 and z < mask.shape[0] and
        y >= 0 and y < mask.shape[1] and
        x >= 0 and x < mask.shape[2]):
        if (mask[z, y, x] == 1):
            mask[z, y, x] = 2
            mask = markregion(mask, (z, y, x+1))
            mask = markregion(mask, (z, y, x-1))    
            mask = markregion(mask, (z, y+1, x))
            mask = markregion(mask, (z, y-1, x))    
            mask = markregion(mask, (z+1, y, x))
            mask = markregion(mask, (z-1, y, x)) 
    return mask

#    from matplotlib import pylab, cm
#    pylab.imshow(image_volume[xyz[2], :, :], interpolation='nearest', origin='lower', cmap=cm.Greys_r, vmin=0, vmax=1000)
#    pylab.plot([roi[0]-0.5, roi[1]-0.5, roi[1]-0.5, roi[0]-0.5, roi[0]-0.5],
#               [roi[2]-0.5, roi[2]-0.5, roi[3]-0.5, roi[3]-0.5, roi[2]-0.5])
#    pylab.scatter(xyz[0]-0.5, xyz[1]-0.5, marker='x')
#    pylab.scatter(xyz_obs[0]-0.5, xyz_obs[1]-0.5, marker='x', color='red')
#    #pylab.scatter(xy_obs[0]-0.5, xy_obs[1]-0.5, marker='x', color='yellow')
#    pylab.show()
#    import numpy
#    total_image = numpy.zeros(dtype=numpy.int32, shape=(roi[3]-roi[2], roi[1]-roi[0]))
#    for z in range(roi[4], roi[5]+1):
#        image = image_volume[z,roi[2]:roi[3],roi[0]:roi[1]]
#        total_image = total_image + image

#    xy = centroid2d(flex.int(total_image), (0, roi[1]-roi[0], 0, roi[3]-roi[2]))
#    print xy[0] + roi[0], xy[1] + roi[2]
#    pylab.imshow(total_image, interpolation='nearest', origin='lower', cmap=cm.Greys_r, vmin=0, vmax=10000)
#    pylab.show()
    
        
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
    usage = "usage: %prog [options] /path/to/GXPARM.XDS /path/to/INTEGRATE.HKL /cbf/search/path/*.cbf"
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
    if len(args) < 3:
        print parser.print_help()
    else:
        create_reflection_mask(args[0], 
                               args[1],
                               args[2], 
                               options.dmin,
                               options.sigma_d,
                               options.sigma_m,
                               options.n_sigma,
                               options.display_frame,
                               options.display_spot)
