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

def visualize_xds_transform_extra(image_volume_original,
                                  image_volume,
                                  xds_grid,
                                  reflections,
                                  spot):
    from matplotlib import pylab, cm
    import numpy
    roi = reflections[spot].region_of_interest
    spot_volume_orig = image_volume_original[roi[4]:roi[5], roi[2]:roi[3], roi[0]:roi[1]]
    spot_volume_bsub = image_volume[roi[4]:roi[5], roi[2]:roi[3], roi[0]:roi[1]]
    spot_grid = xds_grid[reflections[spot].transform_index]

    n_slice_vol_orig = spot_volume_orig.shape[0]
    n_slice_vol = spot_volume_bsub.shape[0]
    n_slice_grid = spot_grid.shape[0]

    blank = numpy.zeros(spot_volume_orig.shape[1:])
    for i in range(0,9):
        ax = pylab.subplot(3, 9, i+1)
        if i < n_slice_vol_orig:
            image = spot_volume_orig[i]
        else:
            image = blank
        plt=pylab.imshow(image, vmin=0, vmax=1000,
            cmap=cm.Greys_r, interpolation='nearest')
        plt.axes.get_xaxis().set_ticks([])
        plt.axes.get_yaxis().set_ticks([])
        ax.set_title("slice: {0}".format(i))

    blank = numpy.zeros(spot_volume_bsub.shape[1:])
    for i in range(0,9):
        ax = pylab.subplot(3, 9, 9+i+1)
        if i < n_slice_vol:
            image = spot_volume_bsub[i]
        else:
            image = blank

        plt=pylab.imshow(image, vmin=0, vmax=1000,
            cmap=cm.Greys_r, interpolation='nearest')
        plt.axes.get_xaxis().set_ticks([])
        plt.axes.get_yaxis().set_ticks([])
        ax.set_title("slice: {0}".format(i))

    for i in range(0,9):
        ax = pylab.subplot(3, 9, 18+i+1)
        image = spot_grid[i]
        plt=pylab.imshow(image, vmin=0, vmax=numpy.max(spot_grid),
            cmap=cm.Greys_r, interpolation='nearest')
        plt.axes.get_xaxis().set_ticks([])
        plt.axes.get_yaxis().set_ticks([])
        ax.set_title("slice: {0}".format(i))
    pylab.show()


def visualize_xds_transform(grid, display_spot):
    """Display the XDS transform for a given spot."""
    from matplotlib import pylab, cm
    import numpy
    spot_grid = grid[display_spot, :, :, :]
#    fig = pylab.figure()
    for i in range(0,9):
        ax = pylab.subplot(3, 3, i+1)
        image = spot_grid[i, :, :]
        plt=pylab.imshow(image, vmin=0, vmax=numpy.max(spot_grid),
            cmap=cm.Greys_r, interpolation='nearest')
        plt.axes.get_xaxis().set_ticks([])
        plt.axes.get_yaxis().set_ticks([])
        ax.set_title("slice: {0}".format(i))
    pylab.show()
#    pylab.savefig("temp/xds_transformed_example.png")
#    for i in range(0,9):
#        #fig = pylab.figure()
#        image = spot_grid[i, :, :]
#        plt=pylab.imshow(image, vmin=0, vmax=numpy.max(spot_grid),
#            cmap=cm.Greys_r)#, interpolation='nearest')
#        plt.axes.get_xaxis().set_ticks([])
#        plt.axes.get_yaxis().set_ticks([])
 #
#        pylab.show()
        #pylab.savefig("temp/xds_transformed_spot_{0:03d}_frame_{1}.png".format(display_spot, i))
        #pylab.clf()

def visualize_xds_transform_3d(grid, display_spot):
    """Display the XDS transform for a given spot."""
    from matplotlib import pylab, cm
    from mpl_toolkits.mplot3d import axes3d
    import numpy
    spot_grid = grid[display_spot, :, :, :]
    gridx = numpy.array([range(spot_grid.shape[2])] * spot_grid.shape[1])
    gridy = numpy.array([range(spot_grid.shape[1])] * spot_grid.shape[2]).transpose()
    minx = 0
    maxx = spot_grid.shape[2]
    miny = 0
    maxy = spot_grid.shape[1]
    minz = numpy.min(spot_grid)
    maxz = numpy.max(spot_grid)
#    for i in range(0,9):
#        fig = pylab.figure()
#        image = spot_grid[i, :, :]
#        ax = fig.add_subplot(111, projection='3d')
#        ax.set_xlim([minx, maxx])
#        ax.set_ylim([miny, maxy])
#        ax.set_zlim([minz, maxz])
#        ax.set_autoscalex_on(False)
#        ax.set_autoscaley_on(False)
#        ax.set_autoscalez_on(False)
#        ax.plot_wireframe(gridx, gridy, image)
#        plt.axes.get_xaxis().set_ticks(map(lambda x: minx + x * (maxx - minx) / 2, range(0, 3)))
#        plt.axes.get_yaxis().set_ticks(map(lambda y: miny + y * (maxy - miny) / 2, range(0, 3)))
#        plt.axes.get_zaxis().set_ticks([])
#        pylab.show()
#        pylab.savefig("temp/xds_transformed_spot_{0}_3d_frame_{1}.png".format(display_spot, i))
#        pylab.clf()

def save_transformed_spot_to_file(grid, display_spot):
    import numpy
    numpy.save('temp/transfomed_spot_{0:03d}'.format(display_spot),
               grid[display_spot, :, :, :])


def get_image_volume_size(list_of_images):

    import os
    import numpy
    from dxtbx.model.detector_helpers_types import detector_helpers_types
    from dxtbx.format.Registry import Registry
    
    for image in list_of_images:
        assert(os.path.exists(image))

    list_of_images.sort()

    format = Registry.find(list_of_images[0])
    
    # verify that these are all the same format i.e. that they are all
    # understood equally well by the format instance.

    format_score = format.understand(list_of_images[0])

    for image in list_of_images:
        assert(format.understand(image) == format_score)

    i = format(list_of_images[0])

    beam = i.get_beam()
    gonio = i.get_goniometer()
    det = i.get_detector()
    scan = i.get_scan()

    # now verify that they share the same detector position, rotation axis
    # and beam properties.

    scans = [scan]

    for image in list_of_images[1:]:
        i = format(image)
        assert(beam == i.get_beam())
        assert(gonio == i.get_goniometer())
        assert(det == i.get_detector())

        scans.append(i.get_scan())

    for s in sorted(scans)[1:]:
        scan += s

    size_z = scan.get_image_range()[1] - scan.get_image_range()[0]
    size_xy = det.get_image_size()

    return (size_z, size_xy[1], size_xy[0])

def get_image_volume(image_paths):

    from iotbx.detectors import ImageFactory
    import numpy
    num_slices, height, width = get_image_volume_size(image_paths)
    
    # Initialise the image volume
    num_slices = len(image_paths)
    volume = numpy.zeros(shape=(num_slices, height, width), dtype=numpy.int32)

    # For each CBF file, read the image and put into the image volume
    for i, filename in enumerate(image_paths):
        image = ImageFactory(filename)
        image.read()
        image_data = image.get_raw_data().as_numpy_array()
        image_data.shape = (height, width)
        volume[i,:,:] = image_data

    # Return image volume
    return volume

def search_for_image_volume(search_path):
    from glob import glob
    filenames = glob(search_path)
    filenames.sort()
    return get_image_volume(filenames)

       
    



def perform_xds_transform(input_filename, cbf_search_path, d_min,
                          sigma_divergence, sigma_mosaicity, n_sigma,
                          display_spot):
    """Read the required data from the file, predict the spots and calculate
       gaussian model parameters."""
    from dials_jmp.spot_prediction import SpotPredictor
    from dials_jmp.integration import ReflectionMaskRoi
    from dials_jmp.integration import ReflectionMask
    from dials_jmp.integration import ReflectionMaskCreator
    from dials_jmp.integration import SubtractBackground
    from dials_jmp.integration import XdsTransformGrid
    from dials_jmp.integration import XdsTransform
    from dials_jmp.io import xdsio, pycbf_extra
    from cctbx import uctbx, sgtbx
    from time import time
    from dials_jmp.array_family import flex
    from dials_jmp.array_family.flex import remove_if_not
    from scitbx import matrix
    from dials_jmp.io.dxtbx_example import factory

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
        #from glob import glob
        #filenames = glob(cbf_search_path)
        #filenames.sort()
        #g, d, b, s, cube = factory(filenames)
        #from iotbx.detectors import ImageFactory
        #image = ImageFactory(self._format % frame_number)
        #image.read()
        #nfast = image.parameters['SIZE2']
        #nslow = image.parameters['SIZE1']
        #frame0 = s.get_image_range()[0]
        #frame1 = s.get_image_range()[1]
        #cube_fast = d.get_image_size()[0]
        #cube_slow = d.get_image_size()[1]
        #print cube_fast, cube_slow
        #image_volume = cube.get(frame0, frame1, 0, cube_slow, 0, cube_fast).as_numpy_array()
        image_volume = search_for_image_volume(cbf_search_path)
        image_volume_original = image_volume.copy()
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

    # Subtract the background for each reflection
    n_reflections = len(region_of_interest)
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

    # Extract arrays from array of reflections
    region_of_interest = flex.tiny6_int(len(reflections))
    image_volume_coords = flex.vec3_double(len(reflections))
    beam_vectors = flex.vec3_double(len(reflections))
    rotation_angles = flex.double(len(reflections))
    for i, r in enumerate(reflections):
        region_of_interest[i] = r.region_of_interest
        image_volume_coords[i] = r.image_coord
        beam_vectors[i] = r.beam_vector
        rotation_angles[i] = r.rotation_angle

    print "Initialising XDS Transform"
    start_time = time()
    n_reflections = len(region_of_interest)
    grid_origin = (4, 4, 4)
    xds_grid = XdsTransformGrid(n_reflections, grid_origin, sigma_divergence,
                                sigma_mosaicity, n_sigma)

    xds_trans = XdsTransform(xds_grid,
                             flex.int(image_volume),
                             reflection_mask_creator.mask,
                             detector, beam, gonio)
    finish_time = time()
    print "Time taken: {0} s".format(finish_time - start_time)

    print "Performing XDS transform on {0} reflections".format(n_reflections)
    start_time = time()
#    xds_trans.calculate(region_of_interest, beam_vectors, rotation_angles)
    xds_trans.calculate(reflections)
    finish_time = time()
    print "Time taken: {0} s".format(finish_time - start_time)


#    index2 = 0
#    for i, r in enumerate(reflections):
#        if (r.mask_index == 29417):
#            index2 = i
#            break
#    r = reflections[index2]
#
#    print r
#
#    display_spot = [r.transform_index]
#
#    roi = r.region_of_interest
#    spot_image = image_volume[roi[4]:roi[5], roi[2]:roi[3], roi[0]:roi[1]]
#    print spot_image

#    r = reflections[3]
#    roi = r.region_of_interest
#    spot_image = image_volume[roi[4]:roi[5], roi[2]:roi[3], roi[0]:roi[1]]
#    print spot_image
#    print r

    if display_spot:
        for spot in display_spot:
            print "Displaying XDS Transform for spot \"{0}\"".format(spot)
            visualize_xds_transform_extra(image_volume_original,
                                          image_volume,
                                          xds_grid.data.as_numpy_array(),
                                          reflections,
                                          spot)
            #visualize_xds_transform(xds_grid.data.as_numpy_array(), spot)
#            visualize_xds_transform_3d(xds_grid.data.as_numpy_array(), spot)
            #save_transformed_spot_to_file(xds_grid.data.as_numpy_array(), spot)
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
