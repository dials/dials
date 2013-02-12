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

def visualize_frame_subtracted_background(image_volume, display_frame, spot_coords):
    """Display the reflection mask for a given frame."""
    from matplotlib import pylab, cm
    #spot_xy = [(x, y) for x, y, z in spot_coords if display_frame <= z < display_frame+1]
    #xcoords, ycoords = zip(*spot_xy)
    fig = pylab.figure()
    plt = pylab.imshow(image_volume[display_frame,:,:], vmin=0, vmax=1000,
        cmap=cm.Greys_r, interpolation="nearest")
    #pylab.scatter(xcoords, ycoords, marker='x')
    plt.axes.get_xaxis().set_ticks([])
    plt.axes.get_yaxis().set_ticks([])
#    pylab.show()
    frame_str = str(display_frame).zfill(3)
    pylab.savefig("temp/background_subtracted_spots_frame_{0}.tiff".format(frame_str))
    pylab.clf()

def visualize_spot_subtracted_background(image_volume, display_spot, image_volume_coords,
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
    pylab.imshow(image_volume, cmap=cm.Greys_r, interpolation="nearest",
        extent=[roi2[0], roi2[1], roi2[2], roi2[3]])
    pylab.plot([roi[0], roi[1], roi[1], roi[0], roi[0]],
               [roi[2], roi[2], roi[3], roi[3], roi[2]])
    pylab.show()

def subtract_reflection_background(input_filename, cbf_search_path, d_min,
                                   sigma_divergence, sigma_mosaicity, n_sigma,
                                   display_frame, display_spot):
    """Read the required data from the file, predict the spots and calculate
       gaussian model parameters."""
    from dials.spot_prediction import SpotPredictor
    from dials.integration import ReflectionMaskRoi
    from dials.integration import ReflectionMask
    from dials.integration import ReflectionMaskCreator
    from dials.integration import SubtractBackground
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

    print image_volume[0,203,236]

    # Extract arrays from array of reflections
    region_of_interest = flex.tiny6_int(len(reflections))
    image_volume_coords = flex.vec3_double(len(reflections))
    for i, r in enumerate(reflections):
        region_of_interest[i] = r.region_of_interest
        image_volume_coords[i] = r.image_coord

    range_x = [roi[1] - roi[0] for roi in region_of_interest]
    range_y = [roi[3] - roi[2] for roi in region_of_interest]
    range_z = [roi[5] - roi[4] for roi in region_of_interest]
    volume = [rx * ry * rz for rx, ry, rz in zip(range_x, range_y, range_z)]
    print "Min/Max ROI X Range: ", min(range_x), max(range_x)
    print "Min/Max ROI Y Range: ", min(range_y), max(range_y)
    print "Min/Max ROI Z Range: ", min(range_z), max(range_z)
    print "Min/Max Roi Volume:  ", min(volume), max(volume)

    index = 28000
    #from dials.integration import BackgroundIntensity
#    import numpy
#    #bi = BackgroundIntensity()
#    roi = region_of_interest[index]
#    spot_data = image_volume[roi[4]:roi[5], roi[2]:roi[3], roi[0]:roi[1]].astype(numpy.float64).copy()
#    #mask_data = reflection_mask_creator.mask[roi[4]:roi[5], roi[2]:roi[3], roi[0]:roi[1]].as_numpy_array()
#    #spot_image = spot_data.copy()[0,:,:]
#    spot_data.shape = -1
#    from matplotlib import pylab, cm
#    spot_data = sorted(spot_data)[::-1]

#    #max_image = numpy.max(spot_image)
#    for i in range(7):
#        fig = pylab.figure()
##        plt = pylab.imshow(spot_image, interpolation='nearest', origin='lower', cmap=cm.Greys_r, vmin=0, vmax=max_image)
#        mean_value = float(numpy.mean(spot_data))
#        #print mean_value
#        #plt.axes.get_xaxis().set_ticks([])
#        #plt.axes.get_yaxis().set_ticks([])
#        print len(spot_data)
#        print spot_data
#        pylab.hist(spot_data)
#        pylab.axvline(mean_value, color='black', linewidth=2,linestyle='--')
#        #pylab.show()
#        filename = "temp/background_calculation_histogram_mean_{0}.png".format(i)
#        print filename
#        pylab.savefig(filename)
#        pylab.clf()
#        #print numpy.mean(spot_data)
#        spot_data = spot_data[1:]
        #ind = numpy.argmax(spot_image)
        #ind = numpy.unravel_index(ind, spot_image.shape)
        #spot_image[ind] = 0
#
#    #print reflections[28000].mask_index
#    #print xyz[28000], xyz[29417]
#    print mask_data
#    print spot_data
#    spot_data = flex.double(spot_data)
#    print "Intensity Value: ", bi.calculate(spot_data)
    print 1/0

    # Subtract the background for each reflection
    print "Subtracting background"
    start_time = time()
    image_volume = flex.int(image_volume)
    subtract_background = SubtractBackground(image_volume, reflection_mask_creator.mask)
    valid_background = flex.bool(len(region_of_interest))
    valid_background = subtract_background.subtract(reflections)
    subtract_background.set_non_reflection_value(0)
    for r, s in zip(reflections, valid_background):
        if (s == False):
            reflection_mask_creator.set_reflection_pixel_value(image_volume, r, 0)
    reflections = remove_if_not(reflections, valid_background)
    image_volume = image_volume.as_numpy_array()
    finish_time = time()
    print "Time taken: {0} s".format(finish_time - start_time)

    index2 = 0
    for i, r in enumerate(reflections):
        if (r.mask_index == 29417):
            index2 = i
            break
    r = reflections[index2]

#    from matplotlib import pylab, cm
#    roi = r.region_of_interest
#    xyz = r.image_coord
#    image_volume = image_volume[xyz[2]-4:xyz[2]+4+1, roi[2]:roi[3], roi[0]:roi[1]]

#    from matplotlib import pylab, cm
#    for z in range(9):
#        fig = pylab.figure()
#        plt = pylab.imshow(image_volume[z,:,:], interpolation='nearest',
#                     origin='lower', cmap=cm.Greys_r,
#                     vmin=0, vmax=numpy.max(image_volume))
#        plt.axes.get_xaxis().set_ticks([])
#        plt.axes.get_yaxis().set_ticks([])
#        pylab.savefig("temp/background_subtracted_spot_intensity_frame_{0}.tiff".format(z))
##        pylab.show()

    #print len(reflections)

    # Extract arrays from array of reflections
    region_of_interest = flex.tiny6_int(len(reflections))
    image_volume_coords = flex.vec3_double(len(reflections))
    for i, r in enumerate(reflections):
        region_of_interest[i] = r.region_of_interest
        image_volume_coords[i] = r.image_coord

#    if display_frame:
#        for frame in display_frame:
#            print "Displaying subtracted background for frame \"{0}\"".format(frame)
#            visualize_frame_subtracted_background(image_volume,
#                frame, image_volume_coords)

#    if display_spot:
#        for spot in display_spot:
#            print "Displaying subtracted background for spot \"{0}\"".format(spot)
#            visualize_spot_subtracted_background(image_volume,
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
        subtract_reflection_background(args[0],
                                       args[1],
                                       options.dmin,
                                       options.sigma_d,
                                       options.sigma_m,
                                       options.n_sigma,
                                       options.display_frame,
                                       options.display_spot)
