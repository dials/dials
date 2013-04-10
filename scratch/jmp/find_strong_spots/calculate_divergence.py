
def calculate_threshold(sweep, trusted_range):
    from scipy.ndimage.measurements import histogram
    from thresholding import maximum_deviation
    import numpy

    threshold_list = []
    for i, flex_image in enumerate(sweep):

        # Get image as numpy array
        image = flex_image.as_numpy_array()

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
        threshold = maximum_deviation(histo)
        threshold_list.append(threshold)

    # Return mean threshold
    return numpy.mean(threshold_list)


def select_strong_pixels(sweep, trusted_range):
    
    import numpy

    # Calculate the threshold
    print "Calculating a threshold."
    threshold = calculate_threshold(sweep, trusted_range)

    # Select only those pixels with counts > threshold
    print "Selecting pixels"
    image = sweep.to_array().as_numpy_array()
    mask = image >= threshold
    return image, mask

def group_pixels(mask):
    from scipy.ndimage.measurements import label, find_objects

    # Label the indices in the mask
    regions, nregions = label(mask)#, structure)

    # Get the bounding box of each object
    objects = find_objects(regions)

    # Return the list of objects
    return objects

def filter_objects(mask, objects, min_pixels):
    import numpy

    # Select objects with more pixels than min_pixels
    new_obj = []
    for obj in objects:
        ind = numpy.where(mask[obj] != 0)[0]
        if len(ind) >= min_pixels:
            new_obj.append(obj)
        else:
            mask[obj] = 0

    return new_obj


def direction_var(values, weights):
    import numpy
    from scitbx import matrix
    weights = numpy.array(weights)
    valx = numpy.array([x for x, y, z in values])
    valy = numpy.array([y for x, y, z in values])
    valz = numpy.array([z for x, y, z in values])

    # Calculate avergae x, y, z
    avrx = numpy.average(valx, weights=weights)
    avry = numpy.average(valy, weights=weights)
    avrz = numpy.average(valz, weights=weights)

    # Calculate mean direction vector
    s1m = matrix.col((avrx, avry, avrz)).normalize()

    # Calculate angles between vectors
    angles = []
    for s in values:
        angles.append(s1m.angle(s))

    # Calculate variance of angles
    angles = numpy.array(angles)
    var = numpy.dot(weights, (angles)**2)/numpy.sum(weights)
    return var

def centroid(image, mask, detector):
    from numpy import zeros, int32, argmax, where, average
    from scitbx.array_family import flex
    from scitbx import matrix

    var = []
    cent = []
    for obj in objects:
        xs = []
        ys = []
        zs = []
        s1s = []
        weights = []
        bbox = [obj[2].start, obj[2].stop,
                obj[1].start, obj[1].stop,
                obj[0].start, obj[0].stop]

        # Calcyulate beam vector for each point
        for i in range(bbox[4], bbox[5]):
            for s in range(bbox[2], bbox[3]):
                for f in range(bbox[0], bbox[1]):
                    if mask[i, s, f]:
                        s1 = matrix.col(detector.get_pixel_lab_coord((f, s)))
                        s1 = s1.normalize()
                        xs.append(f)
                        ys.append(s)
                        zs.append(i)
                        s1s.append(s1)
                        weights.append(image[i, s ,f])

        v = direction_var(s1s, weights)
        var.append(v)
        avrx = average(xs, weights=weights)
        avry = average(ys, weights=weights)
        avrz = average(zs, weights=weights)
        cent.append((avrx, avry, avrz))

    # Return a list of centroids and variances
    return cent, var


def calculate_sigma_beam_divergence(var):
    '''Calculate the beam divergence as the sum of centroid variance of the
    intensity weighted diffracted beam directions.'''
    from math import sqrt

    # Calculate the sum of s^2
    sum_variance = reduce(lambda x, y: x + y, var)

    # Return the beam divergence as the sum / num reflections
    return sqrt(sum_variance / len(var))
    

if __name__ == '__main__':

    import os
    import libtbx.load_env
    from dials.util.nexus import NexusFile 
    from glob import glob 
    from dxtbx.sweep import SweepFactory 
    from math import pi
    
    try:
        dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
        print 'FAIL: dials_regression not configured'
        raise

    # Read the reflections file
    filename = os.path.join(dials_regression, 
        'centroid_test_data', 'reflections.h5')
    handle = NexusFile(filename, 'r')
    
    # Get the reflection list
    print 'Reading reflections.'
    reflections = handle.get_reflections()
    print 'Read {0} reflections.'.format(len(reflections))
    
    # Read images
    template = os.path.join(dials_regression, 
        'centroid_test_data', 'centroid_*.cbf')    
    filenames = glob(template)    
    
    # Load the sweep
    print 'Loading sweep'
    sweep = SweepFactory.sweep(filenames)
    print 'Loaded sweep of {0} images.'.format(len(sweep))
    
    # Select the strong pixels to use in the divergence calculation
    print 'Select the strong pixels from the images.'
    trusted_range = (0, 20000)
    image, mask = select_strong_pixels(sweep, trusted_range)
    
    # Putting pixels into groups
    print 'Putting pixels into groups'
    objects = group_pixels(mask)
    print 'Found {0} objects'.format(len(objects))
    
    print 'Filtering objects'
    min_pixels = 6
    objects = filter_objects(mask, objects, min_pixels)
    print '{0} remaining objects'.format(len(objects))
    
    print 'Calculating centroid and variance.'
    cent, var = centroid(image, mask, sweep.get_detector())
    
    print 'Calculate the beam divergence'
    sigma_d = calculate_sigma_beam_divergence(var)
    print 'Sigma_d = {0} deg'.format(sigma_d * 180.0 / pi)
