#!/usr/bin/env cctbx.python
# cctbx_introduction_1.py
#
# An illustration of how to use cctbx code with imgCIF / cbf files, through
# the pycbf API now included in cctbx. N.B. this should work with the files
# from a Pilatus 300K instrument collected during commissioning on I19 as
# ximg2701_00001.cbf (included).

import math
import sys
import os

# tools to use from CBFlib
import pycbf

# and cctbx
from scitbx import matrix

# This was the known UB matrix, derived from XDS processing.

UB = matrix.sqr([0.0144873, -0.0128813, -0.0002988,
                -0.0128113, -0.0143530, -0.0024004,
                 0.0013736,  0.0019910, -0.0192366])

# Why doesn't Python include a basic nearest integer?!

def nint(a):
    return int(round(a) - 0.5) + (a > 0)

def find_beam_direction(cbf_handle):
    '''Find the beam direction (why is this not simpler in pycbf?)'''

    cbf_handle.find_category('axis')
    cbf_handle.find_column('equipment')
    cbf_handle.find_row('source')

    beam_direction = []

    for j in range(3):
        cbf_handle.find_column('vector[%d]' % (j + 1))
        beam_direction.append(cbf_handle.get_doublevalue())

    B = - matrix.col(beam_direction).normalize()

    return B

def compute_central_rotation_matrix(gonio):
    '''Compute the composite rotation matrix mid way through the acquisition
    of this frame.'''

    x = gonio.rotate_vector(0.5, 1, 0, 0)
    y = gonio.rotate_vector(0.5, 0, 1, 0)
    z = gonio.rotate_vector(0.5, 0, 0, 1)

    R = matrix.rec(x + y + z, (3, 3)).transpose()

    return R

def plot_image(size1, size2, image_values):
    '''Plot an image (values).'''

    from matplotlib import pyplot as plt
    import numpy

    assert(len(image_values) == size1 * size2)

    image = numpy.reshape(numpy.array(image_values, dtype = float),
                          (size1, size2))

    plt.imshow(image)
    plt.savefig('cctbx_introduction_1.png')

    return

def compute_reciprocal_space_distance_map(cbf_image):
    '''Compute a map of distance from the transformed centre of each pixel
    to the nearest reciprocal space node, measured in h, k, l space.'''

    # construct and link a cbf_handle to the image itself.
    cbf_handle = pycbf.cbf_handle_struct()
    cbf_handle.read_file(cbf_image, pycbf.MSG_DIGEST)

    # find the direction if the S0 beam vector (i.e. source -> sample)
    B = find_beam_direction(cbf_handle)

    # find out about the detector
    detector = cbf_handle.construct_detector(0)

    # this returns slow fast slow fast pixels pixels mm mm
    detector_normal = tuple(detector.get_detector_normal())
    distance = detector.get_detector_distance()
    pixel = (detector.get_inferred_pixel_size(1),
             detector.get_inferred_pixel_size(2))

    # find out about the goniometer
    gonio = cbf_handle.construct_goniometer()
    R = compute_central_rotation_matrix(gonio)

    # this method returns slow then fast dimensions i.e. (y, x)

    size = tuple(reversed(cbf_handle.get_image_size(0)))

    wavelength = cbf_handle.get_wavelength()

    O = matrix.col(detector.get_pixel_coordinates(0, 0))
    fast = matrix.col(detector.get_pixel_coordinates(0, 1))
    slow = matrix.col(detector.get_pixel_coordinates(1, 0))

    X = fast - O
    Y = slow - O

    X = X.normalize()
    Y = Y.normalize()
    N = X.cross(Y)

    S0 = (1.0 / wavelength) * B

    # RUBI is (R * UB).inverse()

    RUBI = (R * UB).inverse()

    square_distances = []

    for i in range(size[0]):
        for j in range(size[1]):
            p = matrix.col(detector.get_pixel_coordinates(j, i)).normalize()
            q = (1.0 / wavelength) * p - S0
            _hkl = (RUBI * q).elems
            hkl = map(nint, _hkl)

            d2 = sum([(_hkl[k] - hkl[k]) * (_hkl[k] - hkl[k])
                      for k in range(3)])

            square_distances.append(d2)

    # plot - this needs the full fat cctbx with extra bits

    try:
        plot_image(size[0], size[1], square_distances)
    except:
        print 'Plotting image failed'
 
    # clean up...
 
    detector.__swig_destroy__(detector)
    del(detector)

    gonio.__swig_destroy__(gonio)
    del(gonio)

if __name__ == '__main__':
    compute_reciprocal_space_distance_map(sys.argv[1])
    
