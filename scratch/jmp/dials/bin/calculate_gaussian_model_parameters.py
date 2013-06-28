from __future__ import division
#!/usr/bin/env python

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

def calculate_gaussian_model_parameters(input_filename, cbf_search_path, d_min):
    """Read the required data from the file, predict the spots and calculate
       gaussian model parameters."""
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

    print "Finding strong pixels"
    from dials.array_family.flex import partial_sort_indices
    num_strong_pixels = 100
    flat_image_volume = flex.int(image_volume)
    strong_pixel_indices = partial_sort_indices(flat_image_volume, num_strong_pixels)
    strong_pixel_indices = strong_pixel_indices[0:num_strong_pixels]
    strong_pixel_values = [flat_image_volume[x] for x in strong_pixel_indices]

    print "Assigning indices of nearest reflection"
    from scipy.spatial.kdtree import KDTree
    kd_tree = KDTree(image_volume_coords)
    strong_pixel_xyz = []
    depth, height, width = image_volume.shape
    for i in strong_pixel_indices:
        z = int(i / (width * height))
        zz = z * (width * height)
        y = int((i-zz) / width)
        x = (i-zz) % width
        strong_pixel_xyz.append((x,y,z))
    nearest_neighbours = [kd_tree.query(xyz) for xyz in strong_pixel_xyz]

    print "Sorting the strong pixels by reflection index"
    nearest_neighbours = [int(nn[1]) for nn in nearest_neighbours]
    nearest_neighbours = sorted(enumerate(nearest_neighbours), key=lambda x:x[1])
    strong_pixels = []
    for (i, r) in nearest_neighbours:
        strong_pixels.append((r, strong_pixel_xyz[i], strong_pixel_values[i]))

    print "Find strong reflection ROIs"
    temp = []
    first = True
    for rind, xyz, value in strong_pixels:
        if first:
            last_count = 0
            first = False
        else:
            if (rind != last_rind):
                last_count += 1

        temp.append((last_count, rind, xyz, value))
        last_rind = rind

    strong_pixels = temp
    num_reflections = strong_pixels[-1][0] + 1
    reflection_roi = [[0, 0, 0, 0, False] for i in range(num_reflections)]

    for rcount, rind, (x, y, z), value in strong_pixels:
        roi = reflection_roi[rcount]
        if roi[4] == False:
            roi = [x, x, y, y, True]
        else:
            if x < roi[0]:
                roi[0] = x
            if x > roi[1]:
                roi[1] = x
            if y < roi[2]:
                roi[2] = y
            if y > roi[3]:
                roi[3] = y
        reflection_roi[rcount] = roi

    for r in reflection_roi:
        print r


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
    (options, args) = parser.parse_args()

    # Print help if no arguments specified, otherwise call spot prediction
    if len(args) < 2:
        print parser.print_help()
    else:
        calculate_gaussian_model_parameters(args[0],
                                            args[1],
                                            options.dmin)
