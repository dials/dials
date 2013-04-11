
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
    print "Threshold Value: {0}".format(threshold)

    # Select only those pixels with counts > threshold
    print "Selecting pixels"
    image = sweep.to_array().as_numpy_array()
    mask = image >= threshold
    return image, mask



def find_nearest_neighbour(image, mask, reflections):
    from annlib_ext import AnnAdaptor
    from scitbx.array_family import flex
    from math import sqrt
    import numpy

    # Get the predicted coords
    pred_xyz = []
    for r in reflections:
        x = r.image_coord_px[0]
        y = r.image_coord_px[1]
        z = r.frame_number
        pred_xyz.append((x, y, z))

    # Create the KD Tree
    ann = AnnAdaptor(flex.double(pred_xyz).as_1d(), 3)

    pixel_xyz = []
    ind = numpy.where(mask != 0)
    z = ind[0]
    y = ind[1]
    x = ind[2]

    ann.query(flex.double(zip(x, y, z)).as_1d())

#    for i in xrange(len(ann.nn)):
#        print "Neighbor of {0}, index {1} distance {2}".format(
#        obs_xyz[i], ann.nn[i], sqrt(ann.distances[i]))

    owner = numpy.zeros(shape=mask.shape, dtype=numpy.int32)
    owner[ind] = ann.nn.as_numpy_array()

    return owner

def group_pixels(regions):
    from scipy.ndimage.measurements import find_objects

    # Get the bounding box of each object
    objects = find_objects(regions)

    new_objects = []
    ref_index = []
    for i, obj in enumerate(objects):
        if obj != None:
            new_objects.append(obj)
            ref_index.append(i)

    # Return the list of objects
    return new_objects, ref_index

def filter_objects(mask, objects, ref_index, min_pixels):
    import numpy

    # Select objects with more pixels than min_pixels
    new_obj = []
    new_ind = []
    for obj, rind in zip(objects, ref_index):
        ind = numpy.where(mask[obj] != 0)[0]
        if len(ind) >= min_pixels:
            print len(ind)
            new_obj.append(obj)
            new_ind.append(rind)
        else:
            mask[obj] = 0

    return new_obj, new_ind


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

def centroid(objects, image, mask, detector):
    from numpy import zeros, int32, argmax, where, average
    from scitbx.array_family import flex
    from scitbx import matrix
    import numpy

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
                        s1 = matrix.col(detector.get_pixel_lab_coord((f+0.5, s+0.5)))
                        s1 = s1.normalize()
                        xs.append(f + 0.5)
                        ys.append(s + 0.5)
                        zs.append(i + 0.5)
                        s1s.append(s1)
                        weights.append(image[i, s ,f])

        v = direction_var(s1s, weights)
        var.append(v)
        avrx = average(xs, weights=weights)
        avry = average(ys, weights=weights)
        avrz = average(zs, weights=weights)
        cent.append((avrx, avry, avrz))

    # Return a list of centroids and variances
    return numpy.array(cent), numpy.array(var)


def filter_objects_by_distance(ref_index, reflections, cent, max_distance):

    import numpy
    from math import sqrt

    assert(len(cent) == len(ref_index))

    # Only accept objects closer than max distance
    indices = []
    for i in range(len(ref_index)):
        rx, ry = reflections[ref_index[i]].image_coord_px
        rz = reflections[ref_index[i]].frame_number
        cx, cy, cz = cent[i]
        print rx, ry, rz, cx, cy, cz
        if sqrt((cx - rx)**2 + (cy - ry)**2 + (cz - rz)**2) > max_distance:
            continue
        else:
            indices.append(i)

    return numpy.array(indices)


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

    print 'Find pixel nearest predicted reflection'
    owner = find_nearest_neighbour(image, mask, reflections)

    # Putting pixels into groups
    print 'Putting pixels into groups'
    objects, ref_index = group_pixels(owner)
    print 'Found {0} objects'.format(len(objects))

    print 'Filtering objects'
    min_pixels = 6
    objects, ref_index = filter_objects(mask, objects, ref_index, min_pixels)
    print '{0} remaining objects'.format(len(objects))

    print 'Calculating centroid and variance.'
    cent, var = centroid(objects, image, mask, sweep.get_detector())

    print 'Filter objects by distance from nearest reflection.'
    max_distance = 2
    indices = filter_objects_by_distance(ref_index, reflections, cent, max_distance)
    print '{0} remaining objects'.format(len(indices))
