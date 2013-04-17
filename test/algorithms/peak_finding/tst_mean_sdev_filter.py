
def random_image(size):
    '''Create a random image.'''
    from scitbx.array_family import flex
    image = flex.random_double(size[0] * size[1], 100)
    image.reshape(flex.grid(*size))
    return image

def calculate_mean(image, size):
    '''Calculate the mean filtered image.'''
    from scitbx.array_family import flex
    image_size = image.all()
    mean_image = flex.double(image.accessor())
    num = (2 * size[0] + 1) * (2 * size[1] + 1)
    for j in range(image_size[1]):
        for i in range(image_size[0]):
            i0, i1 = i - size[0], i + size[0] + 1
            j0, j1 = j - size[1], j + size[1] + 1
            if i0 < 0: i0 = 0
            if j0 < 0: j0 = 0
            if i0 > image_size[1]: i0 = image_size[1]
            if j0 > image_size[0]: j0 = image_size[0]

            sub_image = image[j0:j1,i0:i1]
            mean_image[j,i] = flex.sum(sub_image) / num

    return mean_image

def calculate_sdev(image, mean, size):
    '''Calculate the sdev filtered image.'''
    from scitbx.array_family import flex
    from math import sqrt
    import numpy
    image_size = image.all()
    sdev_image = flex.double(image.accessor())
    num = (2 * size[0] + 1) * (2 * size[1] + 1)
    for j in range(image_size[1]):
        for i in range(image_size[0]):
            i0, i1 = i - size[0], i + size[0] + 1
            j0, j1 = j - size[1], j + size[1] + 1
            if i0 < 0: i0 = 0
            if j0 < 0: j0 = 0
            if i0 > image_size[1]: i0 = image_size[1]
            if j0 > image_size[0]: j0 = image_size[0]

            sub_image = image[j0:j1,i0:i1]

            sdev_image[j,i] = numpy.std(sub_image.as_numpy_array())


    return sdev_image

def tst_mean_filter():
    '''Test the mean filter.'''
    from dials.algorithms.peak_finding import mean_filter

    image_size = (256, 256)
    filter_size = (5, 5)

    # Generate a random image
    image = random_image(image_size)

    # Do the filtering
    mean = mean_filter(image, filter_size)

    # Extract the mean from the image
    mean2 = calculate_mean(image, filter_size)

    # Compare values
    for m1, m2 in zip(mean, mean2):
        assert(abs(m1 - m2) < 1e-7)

def tst_mean_sdev_filter():
    '''Test the mean and sdev filter.'''
    from dials.algorithms.peak_finding import mean_sdev_filter

    image_size = (256, 256)
    filter_size = (5, 5)

    # Generate a random image
    image = random_image(image_size)

    # Do the filtering
    mean, sdev = mean_sdev_filter(image, filter_size)

    # Extract the mean from the image
    mean2 = calculate_mean(image, filter_size)
    sdev2 = calculate_sdev(image, mean2, filter_size)

    # Compare values
    for m1, m2 in zip(mean, mean2):
        assert(abs(m1 - m2) < 1e-7)

    # Check the bit inside the edges is ok
    for s1, s2 in zip(sdev[5:-5,5:-5], sdev2[5:-5,5:-5]):
        assert(abs(s1 - s2) < 1e-7)

def run():
    '''Run the tests.'''
    tst_mean_filter()
    tst_mean_sdev_filter()

if __name__ == '__main__':
    run()
