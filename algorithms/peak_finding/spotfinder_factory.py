#!/usr/bin/env python
#
# spot_finder_factory.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division


class FilterRunner(object):
  ''' A class to run multiple filters in succession. '''

  def __init__(self, filters=None):
    ''' Initialise with a list of filters. '''
    if filters is None:
      self.filters=[]
    else:
      self.filters=filters

  def __call__(self, flags, **kwargs):
    ''' Call the filters one by one. '''
    flags = self.check_flags(flags, **kwargs)
    for f in self.filters:
      flags = f(flags, **kwargs)
    return flags

  def check_flags(self, flags, predictions=None, observations=None,
                  shoeboxes=None, **kwargs):
    ''' Check the flags are set, if they're not then create a list
    of Trues equal to the number of items given. '''
    from scitbx.array_family import flex

    # If flags are not set then create a list of Trues
    if flags == None:
      length = 0
      if predictions:
        length = len(predictions)
      if observations:
        if length > 0:
          assert(length == len(observations))
        else:
          length = len(observations)
      if shoeboxes:
        if length > 0:
          assert(length == len(observations))
        else:
          length = len(shoeboxes)

      # Create an array of flags
      flags = flex.bool(length, True)

    # Return the flags
    return flags


class MinPixelsFilter(object):
  ''' Filter the reflections by the number of pixels in the shoeboxes. '''

  def __init__(self, num, code):
    ''' Initialise

    Params:
        num The minimum number of pixels allowed
        code The mask code to use for comparison

    '''
    self.code = code
    self.num = num

  def run(self, flags, observations=None, shoeboxes=None, **kwargs):
    ''' Run the filtering. '''

    # Get the number of mask values matching the code
    count = shoeboxes.count_mask_values(self.code)

    # Return the flags of those > the given number
    return flags.__and__(count >= self.num)

  def __call__(self, flags, **kwargs):
    ''' Call the filter and print information. '''
    from dials.util.command_line import Command
    Command.start('Filtering {0} spots by number of pixels'.format(
        flags.count(True)))
    flags = self.run(flags, **kwargs)
    Command.end('Filtered {0} spots by number of pixels'.format(
        flags.count(True)))
    return flags


class PeakCentroidDistanceFilter(object):

  def __init__(self, maxd):
    ''' Initialise

    Params:
        maxd The maximum distance allowed

    '''
    self.maxd = maxd

  def run(self, flags, observations=None, shoeboxes=None, **kwargs):
    ''' Run the filtering. '''

    # Get the peak locations and the centroids and return the flags of
    # those closer than the min distance
    peak = shoeboxes.peak_coordinates()
    cent = observations.centroids().px_position()
    return flags.__and__((peak - cent).norms() <= self.maxd)

  def __call__(self, flags, **kwargs):
    ''' Call the filter and print information. '''
    from dials.util.command_line import Command
    Command.start('Filtering {0} spots by peak-centroid distance'.format(
        flags.count(True)))
    flags = self.run(flags, **kwargs)
    Command.end('Filtered {0} spots by peak-centroid distance'.format(
        flags.count(True)))
    return flags


class CentroidResolutionFilter(object):

  def __init__(self, d_min, d_max):
    ''' Initialise

    Params:
        dmin The maximum resolution
        dmax The minimum resolution

    '''
    if d_min == None:
      self.d_min = 0.0
    else:
      self.d_min = d_min

    if d_max == None:
      self.d_max = 1000.0
    else:
      self.d_max = d_max

  def run(self, flags, sweep=None, observations=None, **kwargs):
    ''' Run the filtering. '''

    # Get all the observation resolutions
    d = observations.resolution(sweep.get_beam(), sweep.get_detector())

    # Return the flags of those in range
    return (flags.__and__(d >= self.d_min)).__and__(d <= self.d_max)

  def __call__(self, flags, **kwargs):
    ''' Call the filter and print information. '''
    from dials.util.command_line import Command
    Command.start('Filtering {0} spots by resolution'.format(
        flags.count(True)))
    flags = self.run(flags, **kwargs)
    Command.end('Filtered {0} spots by resolution'.format(
        flags.count(True)))
    return flags


class PowderRingFilter(object):

  def __init__(self, crystal_symmetry):
    self.crystal_symmetry = crystal_symmetry

  def run(self, flags, sweep=None, observations=None, **kwargs):
    from cctbx import crystal, sgtbx, uctbx

    from dials.array_family import flex
    from dials.model.data import ReflectionList
    from dxtbx import imageset
    detector = sweep.get_detector()
    beam = sweep.get_beam()

    ms = self.crystal_symmetry.build_miller_set(
      anomalous_flag=False, d_min=detector.get_max_resolution(beam.get_s0()))
    ms = ms.sort(by_value="resolution")

    miller_indices = flex.miller_index()
    two_thetas_obs = flex.double()
    wavelength = beam.get_wavelength()

    for i, centroid in enumerate(observations):
      if not flags[i]: continue
      x, y = centroid.centroid.px_xy
      d_spacing = detector[centroid.panel].get_resolution_at_pixel(
        beam.get_s0(), (x, y))
      for j, d in enumerate(ms.d_spacings().data()):
        if abs(d - d_spacing) < 0.02:
          flags[i] = False
          miller_indices.append(ms.indices()[j])
          two_thetas_obs.append(uctbx.d_star_sq_as_two_theta(
            uctbx.d_as_d_star_sq(d_spacing), wavelength=wavelength, deg=True))

    return flags

  def __call__(self, flags, **kwargs):
    ''' Call the filter and print information. '''
    from dials.util.command_line import Command
    Command.start('Filtering {0} spots by powder rings'.format(
        flags.count(True)))
    flags = self.run(flags, **kwargs)
    Command.end('Filtered {0} spots by powder rings'.format(
        flags.count(True)))
    return flags


class polygon(object):
  def __init__(self, vertices):
    assert len(vertices) > 2
    self.vertices = vertices

  def is_inside(self, x, y):
    # http://en.wikipedia.org/wiki/Point_in_polygon
    # http://en.wikipedia.org/wiki/Even-odd_rule
    poly = self.vertices
    num = len(poly)
    i = 0
    j = num - 1
    inside = False
    for i in range(num):
      if  ((poly[i][1] > y) != (poly[j][1] > y)) and \
          (x < (poly[j][0] - poly[i][0]) * (y - poly[i][1]) / (poly[j][1] - poly[i][1]) + poly[i][0]):
        inside = not inside
      j = i
    return inside


class UntrustedPolygonFilter(object):

  def __init__(self, polygons):
    self.polygons = polygons

  def run(self, flags, sweep=None, observations=None, **kwargs):
    for i, centroid in enumerate(observations):
      if not flags[i]: continue
      x, y = centroid.centroid.px_xy
      for poly in self.polygons:
        if poly.is_inside(x, y):
          flags[i] = False
    return flags

  def __call__(self, flags, **kwargs):
    ''' Call the filter and print information. '''
    from dials.util.command_line import Command
    Command.start('Filtering {0} spots by untrusted polygons'.format(
        flags.count(True)))
    flags = self.run(flags, **kwargs)
    Command.end('Filtered {0} spots by untrusted polygons'.format(
        flags.count(True)))
    return flags


class SpotFinderFactory(object):
  ''' Factory class to create spot finders '''

  @staticmethod
  def from_parameters(params):
    ''' Given a set of parameters, construct the spot finder

    Params:
        params The input parameters

    Returns:
        The spot finder instance

    '''
    from dials.algorithms.peak_finding import SpotFinder

    # Read in the lookup files
    gain_map = SpotFinderFactory.load_image(params.spotfinder.lookup.gain_map)
    dark_map = SpotFinderFactory.load_image(params.spotfinder.lookup.dark_map)
    mask = SpotFinderFactory.load_image(params.spotfinder.lookup.mask)
    params.spotfinder.lookup.gain_map = gain_map
    params.spotfinder.lookup.dark_map = dark_map
    params.spotfinder.lookup.mask = mask
    from dials.framework.registry import Registry
    registry = Registry()
    registry.params().spotfinder = params.spotfinder

    # Configure the algorithm and wrap it up
    find_spots = SpotFinderFactory.configure_algorithm(params)
    filter_spots = SpotFinderFactory.configure_filter(params)
    return SpotFinder(
      find_spots=find_spots,
      filter_spots=filter_spots,
      scan_range=params.spotfinder.scan_range)

  @staticmethod
  def configure_algorithm(params):
    ''' Given a set of parameters, construct the spot finder

    Params:
        params The input parameters

    Returns:
        The spot finder instance

    '''
    from dials.util.command_line import Command
    from dials.algorithms.peak_finding.spot_finder import ExtractSpots

    # Create the threshold strategy
    threshold = SpotFinderFactory.configure_threshold()

    # Setup the spot finder
    return ExtractSpots(threshold_image=threshold,
                        mask=params.spotfinder.lookup.mask)

  @staticmethod
  def configure_threshold():
    ''' Get the threshold strategy'''
    from dials.framework.registry import init_ext
    return init_ext('spotfinder.threshold')

  @staticmethod
  def configure_filter(params):
    ''' Get the filter strategy. '''
    from dials.algorithms import shoebox
    from cctbx import crystal

    # Initialise an empty list of filters
    filters = []

    # Add a min number of pixels filter
    if params.spotfinder.filter.min_spot_size is not None:
      filters.append(MinPixelsFilter(
          params.spotfinder.filter.min_spot_size,
          shoebox.MaskCode.Valid))

    # Add a peak-centroid distance filter
    if params.spotfinder.filter.max_separation is not None:
      filters.append(PeakCentroidDistanceFilter(
        params.spotfinder.filter.max_separation))

    # Add a centroid resolution filter
    if (params.spotfinder.filter.d_min is not None or
        params.spotfinder.filter.d_max is not None):
      filters.append(CentroidResolutionFilter(
          params.spotfinder.filter.d_min,
          params.spotfinder.filter.d_max))

    if params.spotfinder.filter.ice_rings.filter:
      crystal_symmetry = crystal.symmetry(
        unit_cell=params.spotfinder.filter.ice_rings.unit_cell,
        space_group=params.spotfinder.filter.ice_rings.space_group.group())
      filters.append(PowderRingFilter(crystal_symmetry))
    if len(params.spotfinder.filter.untrusted_polygon):
      polygons = []
      for vertices in params.spotfinder.filter.untrusted_polygon:
        if vertices is not None:
          assert len(vertices) % 2 == 0
          vertices = [vertices[i*2:i*2+2] for i in range(len(vertices)//2)]
          polygons.append(polygon(vertices))
      if len(polygons):
        filters.append(UntrustedPolygonFilter(polygons))

    # Return the filter runner with the list of filters
    return FilterRunner(filters)

  @staticmethod
  def load_image(filename):
    ''' Given a filename, load an image

    Params:
        filename The input filename

    Returns:
        The image or None

    '''
    import cPickle as pickle

    # If no filename is set then return None
    if not filename:
      return None

    # Read the image and return the image data
    return pickle.load(open(filename))
