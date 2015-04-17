# -*- coding: utf-8 -*-

"""
Contains the Image class. This class is used to extract all data
necessary to model the experiment from the cbf header. This is done
using dxtbx.
"""

import numpy as np
from numpy.linalg import norm

import dxtbx

import viewer

class Image(object):
    """
    The basic diffraction image class. Used to get useful data from
    a cbf header. Most methods are named so as to be self-explanatory.
    """
    def __init__(self, filename):
        self.image = dxtbx.load(filename)

    def find_beam_direction(self):
        beam = self.image.get_beam()
        direction = beam.get_direction()
        direction /= norm(direction)

        return direction

    def detector_parameters(self):
        """
        Uses dxtbx to get detector parameters. Returns panel
        origins, fast axes, slow axes, panel dimensions, pixel size,
        detector distance, and image resolution.
        """
        detector = self.image.get_detector()

        viewer.num_panels = len(detector)
        d = detector[0].get_distance()

        origin = np.zeros((len(detector), 3))
        fast = np.zeros((len(detector), 3))
        slow = np.zeros((len(detector), 3))
        size = np.zeros((len(detector), 2))
        pix = detector[0].get_pixel_size()
        res = detector[0].get_image_size()

        for i in range(len(detector)):
            origin[i] = detector[i].get_origin()
            fast[i] = detector[i].get_fast_axis()
            slow[i] = detector[i].get_slow_axis()
            size[i] = detector[i].get_image_size_mm()

        return origin, fast, slow, size, pix, d, res

    def get_sample_orientation(self):
        # Currently reading the data using dxtbx, which gives different
        # values to pycbf but may or may not be correct (impossible to
        # know without direct comparison to real beamline images)
        gonio = self.image.get_goniometer()
        fixrot = gonio.get_fixed_rotation()
        i = (fixrot[0], fixrot[3], fixrot[6])
        j = (fixrot[1], fixrot[4], fixrot[7])
        k = (fixrot[2], fixrot[5], fixrot[8])
        # The lines below came from solving linear equations derived
        # from matrix multiplication.
        beta = np.arcsin(-k[0])
        alpha = np.arcsin(k[1] / np.cos(beta))
        gamma = np.arcsin(j[0] / np.cos(beta))
        [alpha, beta, gamma] = np.rad2deg([alpha, beta, gamma])
        axis = gonio.get_rotation_axis()

        return (i, j, k), (alpha, beta, gamma), axis

    def image_to_array(self):
        """
        Converts the diffraction image to a numpy array. Useful
        in other modules when calculating how many pixels lie in
        goniometer shadow.
        """
        raw = self.image.get_raw_data()
        img_array = raw.as_numpy_array()

        return img_array

if __name__ == "__main__":
    image = Image('cbf_image2.cbf')
    a, b, c, d = image.detector_parameters()
    print a, b, c, d
