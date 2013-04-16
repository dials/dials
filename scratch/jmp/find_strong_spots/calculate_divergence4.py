
def calculate_threshold(image, trusted_range):
    from scipy.ndimage.measurements import histogram
    from thresholding import maximum_deviation
    import numpy

    # Cap pixels to within trusted range
    image.shape = -1
    ind = numpy.where(image < trusted_range[0])
    image[ind] = trusted_range[0]
    ind = numpy.where(image > trusted_range[1])
    image[ind] = trusted_range[1]

    # Histogram the pixels
    histo = histogram(image, trusted_range[0], trusted_range[1], trusted_range[1])
    histo = histo / numpy.sum(histo)

    # Calculate the threshold and add to list
    return maximum_deviation(histo)

def select_strong_pixels(sweep, trusted_range):
    from dials.util.command_line import ProgressBar
    import numpy

    # Calculate the threshold
    coordinate = []
    intensity = []
    progress = ProgressBar()
    for i, flex_image in enumerate(sweep):
        image = flex_image.as_numpy_array()
        height, width = image.shape
        threshold = calculate_threshold(image, trusted_range)
        image.shape = -1
        mask = image >= threshold

        ind = numpy.where(mask != 0)[0]
        z = [i] * len(ind)
        y = [int(idx // width) for idx in ind]
        x = [int(idx % width) for idx in ind]
        coords = zip(x, y, z)
        coordinate.extend(coords)
        intensity.extend(list(image[ind]))
        progress.update(100.0 * float(i) / len(sweep))
    progress.finished()

    return coordinate, intensity

def create_groups(pixels, grid_size):

    from scitbx.array_family import flex
    from dials.algorithms.peak_finding import LabelPixels, flex_vec3_int
    label_pixels = LabelPixels(grid_size)
    labels = label_pixels(flex_vec3_int(pixels))
    return labels

if __name__ == '__main__':

    import os
    import libtbx.load_env
    from dials.util.nexus import NexusFile
    from glob import glob
    from dxtbx.sweep import SweepFactory
    from math import pi
    from scitbx import matrix

    try:
        dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
        print 'FAIL: dials_regression not configured'
        raise

    # The XDS values
    xds_sigma_d = 0.060
    xds_sigma_m = 0.154

    # Read the reflections file
    filename = os.path.join(dials_regression,
        'centroid_test_data', 'reflections.h5')
    #filename = os.path.join('/home/upc86896/Data/X1_strong_M1S1_1_', 'reflections.h5')
    handle = NexusFile(filename, 'r')

    # Get the reflection list
    print 'Reading reflections.'
    reflections = handle.get_reflections()
    print 'Read {0} reflections.'.format(len(reflections))

    # Read images
    template = os.path.join(dials_regression,
        'centroid_test_data', 'centroid_*.cbf')
    #template = os.path.join('/home/upc86896/Data/X1_strong_M1S1_1_', 'X1_strong_M1S1_1_*.cbf')
    filenames = glob(template)

    # Load the sweep
    print 'Loading sweep'
    sweep = SweepFactory.sweep(filenames)
    print 'Loaded sweep of {0} images.'.format(len(sweep))

#    # Select the strong pixels to use in the divergence calculation
#    print 'Select the strong pixels from the images.'
#    trusted_range = (0, 20000)
#    coordinate, intensity = select_strong_pixels(sweep, trusted_range)
#    print 'Selected {0} pixels'.format(len(coordinate))
#
#    print 'Create blobs'
#    image_size = sweep.get_detector().get_image_size()
#    grid_size = (image_size[0], image_size[1], len(sweep))
#    labels = create_groups(coordinate, grid_size)
#    print 'Labelled {0} blobs.'.format(max(labels)+1)

    from select_spots import SpotFinder

    print sweep.get_detector()
    sweep.get_detector().set_trusted_range((0, 20000))
    find_spots = SpotFinder()

    find_spots(sweep)
