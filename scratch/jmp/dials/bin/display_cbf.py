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

def display_cbf(cbf_search_path, display_frame):
    """Read the required data from the file, predict the spots and display."""
    from dials.io import xdsio, pycbf_extra
    from time import time

    # Load the image volume from the CBF files and set the number of frames
    print "Searching \"{0}\" for CBF files".format(cbf_search_path)
    image_volume = pycbf_extra.search_for_image_volume(cbf_search_path)

    # If display frame selected then visualize
    for frame in display_frame:
        from matplotlib import pylab, cm
        import numpy
        image = image_volume[frame,:,:]
        print "Displaying CBF file for frame \"{0}\"".format(frame)
        print "Min/Max Image: ", numpy.min(image), numpy.max(image)
        plt = pylab.imshow(image_volume[frame,:,:], vmin=0, vmax=50, cmap=cm.Greys_r)
        plt.axes.get_xaxis().set_ticks([])
        plt.axes.get_yaxis().set_ticks([])
        pylab.show()

if __name__ == '__main__':

    from optparse import OptionParser

    # Specify the command line options
    usage = "usage: %prog [options] /cbf/search/path*.cbf"
    parser = OptionParser(usage)
    parser.add_option('-d', '--display-frame',
                      dest='display_frame',
                      type="string",
                      action="callback",
                      callback=display_frame_callback,
                      help='Select a frame to display with predicted spots')
    (options, args) = parser.parse_args()

    # Print help if no arguments specified, otherwise call spot prediction
    if len(args) == 0:
        print parser.print_help()
    else:display_cbf(args[0],
                     options.display_frame)
        
