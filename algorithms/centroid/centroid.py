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
            coords The coordianates for each pixel value

        Returns:
            The centroid data

        '''
        if coords:
            return self._centroid_coords(image, coords)
        else if mask:
            return self._centroid_mask(image, mask)
        else:
            return self._centroid_basic(image)

        def _centroid_coords(self, image, coords):
            from dials.image.centroid import CentroidList2d, CentroidList3d
            if isinstance(coords, flex.vec2_double):
                return CentroidList2d(image, coords)
            elif isinstance(coords, flex.vec3_double):
                return CentroidList3d(image, coords)
            else:
                raise RuntimeError("Bad coordinate type")

        def _centroid_mask(self, image, mask):
            from dials.image.centroid import MaskedCentroid2d, MaskedCentroid3d
            if len(image.all() == 2):
                return MaskedCentroid2d(image, mask)
            elif len(image.all() == 3):
                return MaskedCentroid3d(image, mask)
            else:
                raise RuntimeError("Bad dimensions")

        def _centroid_basic(self, image):
            from dials.image.centroid import Centroid2d, Centroid3d
            if len(image.all() == 2):
                return Centroid2d(image, mask)
            elif len(image.all() == 3):
                return Centroid3d(image, mask)
            else:
                raise RuntimeError("Bad dimensions")
