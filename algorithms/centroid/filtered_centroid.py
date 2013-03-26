#
# filtered_centroid.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

from dials.interfaces.centroid.centroid_interface_prototype import \
    centroid_interface_prototype as centroid_interface

from dials.interfaces.centroid.centroid_interface_prototype import \
    CentroidException

def select_contigious(index, shape, conn):
    """ Filter the indices to only contain contigious pixels

    Params:
        index The input array of indices
        shape The shape of the mask
        conn The connectivity (4 or 8)

    Returns:
        The filtered indices

    """
    from scipy.ndimage.measurements import label, histogram
    from numpy import zeros, int32, argmax, where
    from operator import mul
    from scitbx.array_family import flex

    # Set the label structure
#    if conn == 8:
#        structure = [[1, 1, 1], [1, 1, 1], [1, 1, 1]]
#    else:
#        structure = [[0, 1, 0], [1, 1, 1], [0, 1, 0]]

    # Create a mask of the roi size
    mask = zeros(reduce(mul, shape), dtype=int32)

    # Check some foreground pixels are given
    if len(index) == 0:
        raise CentroidException("No foreground pixels specified")

    # Set all the points at the indices to 1
    for i in index:
        mask[i] = 1

    # Reshape the mask
    mask.shape = shape

    # Label the indices in the mask
    regions, nregions = label(mask)#, structure)

    # Histogram the regions to get the largest
    histo = histogram(regions, 0, nregions+1, nregions+1)

    # Find largest region
    max_index = argmax(histo[1:]) + 1
    
    # Return only those region indices that are equal to max_index
    regions.shape = (-1)
    return flex.size_t(where(regions == max_index)[0])


def select_foreground_pixels(pixel_data, min_pixels=10, n_sigma=-1, conn=4):
    """ Select the pixels contributing to the spot foreground

    Fixed select all foreground pixels that leave the remaining background
    pixels in an approximate normal distribution. Then filter indices to
    get only those pixels that are contigious

    Params:
        pixel_data The pixel data
        min_pixels The minimum number of pixels contributing to background
        n_sigma The number of standard deviations to check for
        conn The connectivity (4 or 8)

    Returns:
        A list of indices deemed to contribute to spot foreground

    """
    from dials.algorithms.integration import foreground_pixels
    from scitbx.array_family import flex

    # Get the 1D array
    data = flex.double(len(pixel_data))
    for i in range(len(pixel_data)):
        data[i] = pixel_data[i]

    # Get a list of foreground pixel indices
    index = foreground_pixels(data, min_pixels, n_sigma)

    # Select and return only those pixels that are contigious
    return select_contigious(index, pixel_data.all(), conn)


class FilteredCentroid(centroid_interface):
    """ Calculate the centroid filtered by foreground indices """
    
    def __init__(self, reflections, min_pixels=10, n_sigma=3, conn=4):
        """ Initialise the class """

        # Save some parameters needed for selection of pixels
        self.min_pixels = min_pixels
        self.n_sigma = n_sigma
        self.conn = conn

        # Init the interface base class
        centroid_interface.__init__(self, reflections)

    def compute_shoebox_centroid(self, shoebox):
        """ Compute the centroid """
        from scitbx.array_family import flex
        import math

        # Get the size of the shoebox
        f_size, r_size, c_size = shoebox.all()

        # Check shoebox values are valid
        for i in shoebox:
            if i < 0:
                raise CentroidException, 'negative pixels in cube'

        # Get a list of pixels
        pixel_list = []
        for f in range(f_size):
            for r in range(r_size):
                for c in range(c_size):
                    pixel_list.append(
                        (f + 0.5, r + 0.5, c + 0.5, shoebox[f, r, c]))

        # Select the indices to use in centroid calculation
        pixel_index = select_foreground_pixels(shoebox, self.min_pixels,
                                               self.n_sigma, self.conn)

        # compute averages of positions
        f_tot, r_tot, c_tot, d_tot = 0.0, 0.0, 0.0, 0.0
        for i in pixel_index:
            f, r, c, d = pixel_list[i]
            f_tot += d * f
            r_tot += d * r
            c_tot += d * c
            d_tot += d

        # Ensure d_tot is not zero
        if not d_tot:
            raise CentroidException("Invalid value for total intensity")

        # Calculate the centroid position
        _f, _r, _c = f_tot / d_tot, r_tot / d_tot, c_tot / d_tot

        # now compute the variance
        f_tot, r_tot, c_tot, d_tot = 0.0, 0.0, 0.0, 0.0
        for i in pixel_index:
            f, r, c, d = pixel_list[i]
            f_tot += d * (f - _f) ** 2
            r_tot += d * (r - _r) ** 2
            c_tot += d * (c - _c) ** 2
            d_tot += d

        # Compute standard deviation of the centroid position
        _sf = math.sqrt(f_tot / d_tot)
        _sr = math.sqrt(r_tot / d_tot)
        _sc = math.sqrt(c_tot / d_tot)

        # Return the position and standard deviation
        return _f, _r, _c, _sf, _sr, _sc
