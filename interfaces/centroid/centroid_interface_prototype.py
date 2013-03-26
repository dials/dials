#
# centroid_interface_prototype.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: Graeme Winter
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
import exceptions


class CentroidException(exceptions.Exception):
    pass


class centroid_interface_prototype(object):
    '''Centroid calculation interface: all of the stuff in here is in pixels.'''

    def __init__(self, reflections):
        '''Initialise the algorithm with the list of reflections.'''

        from dials.model.data import ReflectionList

        # Save the reflection list
        self._reflections = ReflectionList()

        # Calculate the centroids
        self._compute_centroids(reflections)

    def get_reflections(self):
        '''Return the list of reflections'''
        return self._reflections

    def compute_shoebox_centroid(self, shoebox):
        '''Overload me: shoebox has the form of the subset of data which
        should contain all pixels identified within the reflection bounding
        box.'''

        raise RuntimeError, 'overload me'

    def _compute_centroids(self, reflections):
        ''' Compute the centroids. For each reflection, call the overloaded
        compute_centroid_from_bbox method.'''

        # Loop through all the refflections and try to calculate the
        # centroid. If a CentroidException is encountered, then pass.
        for i, ref in enumerate(reflections):

            try:

                # Compute the centroid
                f, r, c, sf, sr, sc = self.compute_shoebox_centroid(ref.shoebox)

                # Add the bounding box offset to the centroid position
                f += ref.bounding_box[4]
                r += ref.bounding_box[2]
                c += ref.bounding_box[0]

                # Add centroid data to reflection
                ref.centroid_position = (c, r, f)
                ref.centroid_variance = (sc, sr, sf)

                # Copy reflection back into array
                self._reflections.append(ref)

            except CentroidException, e:
                continue

        return
