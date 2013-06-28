from __future__ import division


def calculate_sigma_beam_divergence(reflections):
    '''Calculate the beam divergence as the sum of centroid variance of the
    intensity weighted diffracted beam directions.'''
    from math import sqrt

    # Calculate the sum of s^2
    beam_vector_variance = [r.beam_vector_variance for r in reflections]
    sum_variance = reduce(lambda x, y: x + y, beam_vector_variance)

    # Return the beam divergence as the sum / num reflections
    return sqrt(sum_variance / len(reflections))


def round_up(x, base=1):
    return int(base * round(float(x)/base))


def wvar(values, weights):
    import numpy
    from math import sqrt
    weights = numpy.array(weights)
    valx = numpy.array([x for x, y, z in values])
    valy = numpy.array([y for x, y, z in values])
    valz = numpy.array([z for x, y, z in values])
    avrx = numpy.average(valx, weights=weights)
    avry = numpy.average(valy, weights=weights)
    avrz = numpy.average(valz, weights=weights)
    varx = numpy.dot(weights, (valx-avrx)**2)/numpy.sum(weights)
    vary = numpy.dot(weights, (valy-avry)**2)/numpy.sum(weights)
    varz = numpy.dot(weights, (valz-avrz)**2)/numpy.sum(weights)

    return sqrt(varx**2 + vary**2 + varz**2)

def wcovar(values, weights):
    import numpy
    weights = numpy.array(weights)
    valx = numpy.array([x for x, y, z in values])
    valy = numpy.array([y for x, y, z in values])
    valz = numpy.array([z for x, y, z in values])
    avrx = numpy.average(valx, weights=weights)
    avry = numpy.average(valy, weights=weights)
    avrz = numpy.average(valz, weights=weights)

    #c = 1.0 / (numpy.sum(weights)**2 - numpy.sum(weights**2))
#    c = numpy.sum(weights) / (numpy.sum(weights)**2 - numpy.sum(weights**2))
    c = 1.0 / numpy.sum(weights)
    covxx = c * numpy.sum(weights * (valx - avrx) * (valx - avrx))
    covyy = c * numpy.sum(weights * (valy - avry) * (valy - avry))
    covzz = c * numpy.sum(weights * (valz - avrz) * (valz - avrz))
    covxy = c * numpy.sum(weights * (valx - avrx) * (valy - avry))
    covxz = c * numpy.sum(weights * (valx - avrx) * (valz - avrz))
    covyz = c * numpy.sum(weights * (valy - avry) * (valz - avrz))
    #print covxx
#    return covxx + covyy + covzz + 2*abs(covxy) + 2*abs(covxz) + 2*abs(covyz)
    return covxx + covyy + covzz + 2*covxy + 2*covxz + 2*covyz

#http://smb.slac.stanford.edu/facilities/software/xds/
#
#BEAM_DIVERGENCE = arctan(spot diameter/DETECTOR_DISTANCE))

def avar(values, weights):
    import numpy
    from scitbx import matrix
    weights = numpy.array(weights)
    valx = numpy.array([x for x, y, z in values])
    valy = numpy.array([y for x, y, z in values])
    valz = numpy.array([z for x, y, z in values])
    avrx = numpy.average(valx, weights=weights)
    avry = numpy.average(valy, weights=weights)
    avrz = numpy.average(valz, weights=weights)

    s1m = matrix.col((avrx, avry, avrz)).normalize()
    angles = []
    for s in values:
        angles.append(s1m.angle(s))


    angles = numpy.array(angles)
    var = numpy.dot(weights, (angles)**2)/numpy.sum(weights)
    return var


def centroid(image, mask, detector):
    from scipy.ndimage.measurements import label, histogram, find_objects
    from numpy import zeros, int32, argmax, where
    from operator import mul
    from scitbx.array_family import flex
    from scitbx import matrix


    num_on = numpy.sum(mask)

    # Label the indices in the mask
    regions, nregions = label(mask)#, structure)

    # Get the bounding box of each object
    objects = find_objects(regions)




    var = []
    count = 0
    nel = []
    for obj in objects:
        s1s = []
        weights = []
        bbox = [obj[2].start, obj[2].stop,
                obj[1].start, obj[1].stop,
                obj[0].start, obj[0].stop]

        for i in range(bbox[4], bbox[5]):
            for s in range(bbox[2], bbox[3]):
                for f in range(bbox[0], bbox[1]):
                    if mask[i,s,f]:
                        s1 = matrix.col(detector.get_pixel_lab_coord((f, s)))
                        s1 = s1.normalize()
                        s1s.append(s1)
                        weights.append(image[i, s ,f])

        if len(weights) <= 6:
            continue

        nel = len(weights)
        count += 1
        #print count
        #temp = wvar(s1s, weights)
        #temp = wcovar(s1s, weights)
        temp = avar(s1s, weights)
        #print temp
        var.append(temp)

    print "Mean len: ", numpy.mean(nel)

    from math import sqrt, pi

    # Calculate the sum of s^2
    sum_variance = reduce(lambda x, y: x + y, var)

    # Return the beam divergence as the sum / num reflections
    sigma_d = sum_variance / count
    sigma_d = sqrt(sigma_d)
    sigma_d = sigma_d * 180 / pi
    print sigma_d

if __name__ == '__main__':
    from dxtbx.sweep import SweepFactory
    from glob import glob
    from scipy.ndimage.measurements import histogram
    import numpy
    from scitbx.array_family import flex
    from thresholding import maximum_deviation
    from thresholding import percentage_threshold
    from thresholding import normal_threshold

    #filenames = glob("/home/upc86896/Data/X1_strong_M1S1_1_/X1_strong_M1S1_1_*.cbf")
    filenames = glob("/home/upc86896/Projects/cctbx/sources/dials_regression/centroid_test_data/centroid_000*.cbf")
    sweep = SweepFactory.sweep(filenames)
    trusted_range = (0, 20000)

    histo = []
    threshold_list = []
    start = 0
    stop = 3
    for i, flex_image in enumerate(sweep[start:stop]):

        image = flex_image.as_numpy_array()

        image.shape = -1
        ind = numpy.where(image < trusted_range[0])
        image[ind] = trusted_range[0]
        ind = numpy.where(image > trusted_range[1])
        image[ind] = trusted_range[1]

        mean = numpy.mean(image)
        sdev = numpy.std(image)
        var = sdev **2

        temp = histogram(image, trusted_range[0], trusted_range[1], trusted_range[1])
        temp = temp / numpy.sum(temp)
        histo.append(temp)

        #threshold = bhthreshold(temp)
        threshold = maximum_deviation(temp)
        #threshold = percentage_threshold(temp)
#        threshold = normal_threshold(temp)
        threshold_list.append(threshold)


        print "{0} - Mean: {1}, Sdev: {2}, Var: {3}, Thres: {4}".format(
            i, mean, sdev, var, threshold)

    print "Threshold - Mean: {0}, Sdev: {1}".format(
        numpy.mean(threshold_list), numpy.std(threshold_list))

    mean_thresh = numpy.mean(threshold_list)

    image = sweep[start:stop].to_array().as_numpy_array()
    mask = image > mean_thresh
    centroid(image, mask, sweep.get_detector())
#    for i in range(1):

        #from matplotlib import pylab, cm
        #pylab.plot(histo[i])
        #pylab.show()

        #pylab.imshow(sweep[i].as_numpy_array(), cmap=cm.Greys_r, vmax=mean_thresh)
        #pylab.show()

        #mask = sweep[i].as_numpy_array() > mean_thresh

        #centroid(sweep[i].as_numpy_array(), mask, sweep.get_detector())

        #pylab.imshow(mask, cmap=cm.Greys_r)
        #pylab.show()


#    for h in histo:
#        pylab.plot(h)
#    pylab.show()

#    n = 0
#    mean = 0
#    M2 = 0
#    histo = numpy.zeros(shape=(trusted_range[1]+1), dtype=numpy.int32)
#    for i, image in enumerate(sweep[0:1]):
#        np_image = image.as_numpy_array()
#        np_image.shape = -1
#        ind = numpy.where(np_image < trusted_range[0])
#        np_image[ind] = trusted_range[0]
#        ind = numpy.where(np_image > trusted_range[1])
#        np_image[ind] = trusted_range[1]
#        print i
#        for p in np_image:
#            histo[p] += 1
#            #n = n + 1
#            #delta = p - mean
#            #mean = mean + delta / n
#            #M2 = M2 + delta * (p - mean)

#        print i

#    variance = M2/(n - 1)

#    histo = histo / numpy.sum(histo)

#    histo = histogram(image, 0, 20000, 20000)
#    histo = histo / numpy.sum(histo)
##    print "Mean: ", mean
##    print "Var: ", variance
##
##    from math import log
##    entropy = 0
##    for h in histo:
##        if h > 0:
##            entropy += h * log(h)
##    print -entropy
##
#    from matplotlib import pylab
#    #pylab.plot(histo)
#    #pylab.show()
##
#    cumulative = histo
#    for i in range(1, len(histo)):
#        cumulative[i] += cumulative[i-1]
##
#    pylab.plot(cumulative)
#    pylab.show()
