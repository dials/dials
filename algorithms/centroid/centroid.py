#
# centroid.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

class Centroid(object):
    '''Class to centroid data'''

    def __init__(self, pixels, mask=None, coords=None):
        '''Centroid the data.

        If coords is set then the centroid of the list of coordinates is done.
        Otherwise, if mask is set the the centroid of the masked image is done.
        Otherwise, a simple centroid of the image is done.

        Params:
            pixels The image (or list of pixel values) to centroid
            mask The mask of the image (optional)
            coords The coordianates for each pixel value (optional)

        Returns:
            The centroid data

        '''
        # Choose the algorithm
        if coords:
            result = self._centroid_coords(image, coords)
        elif mask:
            result = self._centroid_mask(image, mask)
        else:
            result = self._centroid_basic(image)

        # Extract the result quantities
        self.position = result.position
        self.variance = result.variance
        self.sq_width = result.sq_width

        def _centroid_coords(self, image, coords):
            '''Do the centroid on a list of coordinates.'''
            from dials.image.centroid import centroid_coord_list
            return centroid_coord_list(image, coords)

        def _centroid_mask(self, image, mask):
            '''Do the centroid on an image and mask.'''
            from dials.image.centroid import MaskedCentroid2d, MaskedCentroid3d
            if len(image.all() == 2):
                return MaskedCentroid2d(image, mask)
            elif len(image.all() == 3):
                return MaskedCentroid3d(image, mask)
            else:
                raise RuntimeError("Bad dimensions")

        def _centroid_basic(self, image):
            '''Do the centroid on an image.'''
            from dials.image.centroid import Centroid2d, Centroid3d
            if len(image.all() == 2):
                return Centroid2d(image, mask)
            elif len(image.all() == 3):
                return Centroid3d(image, mask)
            else:
                raise RuntimeError("Bad dimensions")


def centroid_px_to_mm(detector, scan, position, variance, sq_width):
    '''Convenience function to calculate centroid in mm/rad from px'''
    from operator import mul

    # Get the pixel to millimeter function
    pixel_size = detector.get_pixel_size()
    oscillation = scan.get_oscillation(deg=False)
    scale = pixel_size + (oscillation[1],)
    scale2 = map(mul, scale, scale)

    # Convert Pixel coordinate into mm/rad
    x, y, z = position
    xy_mm = detector.pixel_to_millimeter((x, y))
    z_rad = scan.get_angle_from_array_index(z, deg=False)

    # Set the position, variance and squared width in mm/rad
    # N.B assuming locally flat pixel to millimeter transform
    # for variance calculation.
    position_mm = xy_mm + (z_rad,)
    variance_mm = map(mul, variance, scale2)
    sq_width_mm = map(mul, sq_width, scale2)

    # Return the stuff in mm/rad
    return position_mm, variance_mm, sq_width_mm


class CentroidPxAndMM(object):
    '''Class to calculate centroid in pixels and millimeters'''

    def __init__(self, detector, scan, pixels, mask=None, coords=None):
        '''Calculate the centroid in pixels and millimeters.'''

        # Calculate the centroid and extract the values
        centroid = Centroid(pixels, mask, coords)
        self.position_px = centroid.position
        self.variance_px = centroid.variance
        self.sq_width_px = centroid.sq_width

        # Convert Pixel coordinate into mm/rad
        result = centroid_px_to_mm(detector, scan, position, variance, sq_width)
        self.position_mm, self.variance_mm, self.sq_width_mm = result
