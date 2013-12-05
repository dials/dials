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

    # Read in the lookup files
    gain_map = SpotFinderFactory.load_image(params.lookup.gain_map)
    dark_map = SpotFinderFactory.load_image(params.lookup.dark_map)
    mask = SpotFinderFactory.load_image(params.lookup.mask)

    # Create the threshold strategy
    threshold = SpotFinderFactory.configure_threshold(
      params, gain_map, mask)

    # Setup the spot finder
    return ExtractSpots(threshold_image=threshold)

  @staticmethod
  def configure_threshold(params, gain_map, mask):
    ''' Get the threshold strategy'''
    from dials.algorithms.peak_finding.threshold \
        import UnimodalThresholdStrategy, XDSThresholdStrategy

    # Chose the strategy
    if params.spotfinder.threshold.algorithm == 'xds':
      return XDSThresholdStrategy(
          kernel_size=params.spotfinder.threshold.kernel_size,
          gain=gain_map,
          mask=mask,
          n_sigma_b=params.spotfinder.threshold.sigma_background,
          n_sigma_s=params.spotfinder.threshold.sigma_strong,
          min_count=params.spotfinder.threshold.min_local)

    elif params.spotfinder.threshold.algorithm == 'unimodal':
      return UnimodalThresholdStrategy()

    else:
      raise RuntimeError('Unknown threshold strategy')

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
    from dials.util import image

    # If no filename is set then return None
    if not filename:
      return None

    # Read the image and return the image data
    handle = image.reader()
    handle.read_file(filename)
    return handle.get_data()
