#
# spot_matcher.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division


class SpotMatcher(object):
  '''Match the observed with predicted spots.'''

  def __init__(self, max_separation=2):
    '''
    Setup the algorithm

    :param max_separation: Max pixel dist between predicted and observed spot

    '''
    # Set the algorithm parameters
    self._max_separation = max_separation

  def __call__(self, observed, predicted):
    '''i
    Match the observed reflections with the predicted.

    :param observed: The list of observed reflections.
    :param predicted: The list of predicted reflections.

    :returns: The list of matched reflections

    '''
    from dials.array_family import flex

    # Find the nearest neighbours and distances
    # Command.start('Finding nearest neighbours')
    nn, dist = self._find_nearest_neighbours(observed, predicted)
    # Command.end('Found nearest neighbours')

    # Filter the matches by distance
    # Command.start('Filtering matches by distance')
    index = self._filter_by_distance(nn, dist)
    # Command.end('Filtered {0} matches by distance'.format(len(index)))

    # Filter out duplicates to just leave the closest pairs
    # Command.start('Removing duplicate matches')
    len_index = len(index)
    index = self._filter_duplicates(index, nn, dist)
    len_diff = len_index - len(index)
    # Command.end('Removed {0} duplicate match(es)'.format(len_diff))

    # Copy all of the reflection data for the matched reflections
    return flex.size_t(index), flex.size_t([nn[i] for i in index])

  def _find_nearest_neighbours(self, observed, predicted):
    '''
    Find the nearest predicted spot to the observed spot.

    :param observed: The observed reflections
    :param predicted: The predicted reflections

    :returns: (nearest neighbours, distance)

    '''
    from scitbx.array_family import flex

    # Get the predicted coordinates
    predicted_panel = predicted['panel']
    predicted_xyz = predicted['xyzcal.px']
    observed_panel = observed['panel']
    observed_xyz = observed['xyzobs.px.value']

    # Get the number of panels
    max_panel1 = flex.max(predicted_panel)
    max_panel2 = flex.max(observed_panel)
    max_panel = max([max_panel1, max_panel2])

    nn_all = flex.size_t()
    dd_all = flex.double()
    for panel in range(max_panel+1):
      pind = predicted_panel == panel
      oind = observed_panel == panel
      pxyz = predicted_xyz.select(pind)
      oxyz = observed_xyz.select(oind)
      try:
        nn, d = self._find_nearest_neighbours_single(oxyz, pxyz)
        indices = flex.size_t(range(len(pind))).select(pind)
        indices = indices.select(flex.size_t(list(nn)))
        nn_all.extend(indices)
        dd_all.extend(d)
      except Exception:
        warn("Unable to match spots on panel %d" % panel)
    return nn_all, dd_all

  def _find_nearest_neighbours_single(self, oxyz, pxyz):
    '''
    Find the nearest predicted spot to the observed spot.

    :param observed: The observed reflections
    :param predicted: The predicted reflections

    :returns: (nearest neighbours, distance)

    '''
    from annlib_ext import AnnAdaptor
    from scitbx.array_family import flex

    # Create the KD Tree
    ann = AnnAdaptor(pxyz.as_double().as_1d(), 3)

    # Query to find all the nearest neighbours
    ann.query(oxyz.as_double().as_1d())

    # Return the nearest neighbours and distances
    return ann.nn, flex.sqrt(ann.distances)

  def _filter_by_distance(self, nn, dist):
    '''
    Filter the matches by distance.

    :param nn: The nearest neighbour list
    :param dist: The distances

    :returns: A reduced list of nearest neighbours

    '''
    from scitbx.array_family import flex
    index = range(len(nn))
    return flex.int([i for i in index if dist[i] <= self._max_separation])

  def _filter_duplicates(self, index, nn, dist):
    '''
    Filter the matches to remove duplicates

    :param index: The indices of valid spots
    :param nn: The nearest neighbour indices
    :param dist: The distances

    :returns: A reduced list of nearest neighbours

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
