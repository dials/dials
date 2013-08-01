#
# spot_matcher.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

class SpotMatcher(object):
    '''Match the observed with predicted spots.'''

    def __init__(self, max_separation=2):
      '''Setup the algorithm

      Params:
          max_separation Max pixel dist between predicted and observed spot

      '''
      # Set the algorithm parameters
      self._max_separation = max_separation

    def __call__(self, observed, predicted):
        '''Match the observed reflections with the predicted.

        Params:
            observed The list of observed reflections.
            predicted The list of predicted reflections.

        Returns:
            The list of matched reflections

        '''
        from dials.model.data import ReflectionList
        from dials.util.command_line import Command

        # Find the nearest neighbours and distances
        Command.start('Finding nearest neighbours')
        nn, dist = self._find_nearest_neighbours(observed, predicted)
        Command.end('Found nearest neighbours')

        # Filter the matches by distance
        Command.start('Filtering matches by distance')
        index = self._filter_by_distance(nn, dist)
        Command.end('Filtered {0} matches by distance'.format(len(index)))

        # Filter out duplicates to just leave the closest pairs
        Command.start('Removing duplicate matches')
        len_index = len(index)
        index = self._filter_duplicates(index, nn, dist)
        len_diff = len_index - len(index)
        Command.end('Removed {0} duplicate match(es)'.format(len_diff))

        # Copy all of the reflection data for the matched reflections
        reflections = ReflectionList()
        for i in index:
            o = observed[i]
            p = predicted[nn[i]]
            o.miller_index = p.miller_index
            o.rotation_angle = p.rotation_angle
            o.beam_vector = p.beam_vector
            o.image_coord_px = p.image_coord_px
            o.image_coord_mm = p.image_coord_mm
            o.panel_number = p.panel_number
            o.frame_number = p.frame_number
            reflections.append(o)

        # Return the list of reflections
        return reflections

    def _find_nearest_neighbours(self, observed, predicted):
        '''Find the nearest predicted spot to the observed spot.

        Params:
            observed The observed reflections
            predicted The predicted reflections

        Returns:
            (nearest neighbours, distance)

        '''
        from annlib_ext import AnnAdaptor
        from scitbx.array_family import flex
        from math import sqrt

        # Get the predicted coordinates
        predicted_xyz = []
        for r in predicted:
            x, y = r.image_coord_px
            z = r.frame_number
            predicted_xyz.append((x, y, z))

        observed_xyz = [r.centroid_position for r in observed]

        # Create the KD Tree
        ann = AnnAdaptor(flex.double(predicted_xyz).as_1d(), 3)

        # Query to find all the nearest neighbours
        ann.query(flex.double(observed_xyz).as_1d())

        # Return the nearest neighbours and distances
        return ann.nn, flex.sqrt(ann.distances)

    def _filter_by_distance(self, nn, dist):
        '''Filter the matches by distance.

        Params:
            nn The nearest neighbour list
            dist The distances

        Returns:
            A reduced list of nearest neighbours

        '''
        from scitbx.array_family import flex
        index = range(len(nn))
        return flex.int([i for i in index if dist[i] <= self._max_separation])

    def _filter_duplicates(self, index, nn, dist):
        ''' Filter the matches to remove duplicates

        Params:
            index The indices of valid spots
            nn The nearest neighbour indices
            dist The distances

        Returns:
            A reduced list of nearest neighbours

        '''
        seen = {}
        for i in index:
            p = nn[i]
            if p in seen:
                j = seen[p]
                if dist[i] < dist[j]:
                    seen[p] = i
            else:
                seen[p] = i;

        index = []
        for k, v in seen.iteritems():
            index.append(v)
        return index
