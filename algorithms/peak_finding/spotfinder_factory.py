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

class SpotFinder(object):
    ''' A class to do spot finding and filtering. '''

    def __init__(self, find_spots=None, filter_spots=None, scan_range=None):
        ''' Initialise the class. '''

        # Set the spot finding and filter functions
        assert(find_spots != None and filter_spots != None)
        self.find_spots = find_spots
        self.filter_spots = filter_spots

        # Set the scan range
        self.scan_range = scan_range

    def __call__(self, sweep):
        ''' Do the spot finding '''
        from dials.model.data import ReflectionList
        from dials.array_family import flex
        from dials.util.command_line import Command

        # Get list of scan ranges
        if not self.scan_range:
            scan_range = [(0, len(sweep))]
        else:
            scan_range = self.scan_range

        # Get spots from bits of scan
        spots_all = []
        for scan in scan_range:
            j0, j1 = scan
            assert(j1 > j0 and j0 >= 0 and j1 <= len(sweep))
            print '\nFinding spots in image {0} to {1}...'.format(j0, j1)

            # Find the spots
            spots = self.find_spots(sweep[j0:j1])

            # Add the spots to the list
            spots_all.extend(spots)

        # Get the list of shoeboxes
        shoeboxes = flex.shoebox(spots_all)

        # Calculate the spot centroids
        Command.start('Calculating {0} spot centroids'.format(len(shoeboxes)))
        centroid = shoeboxes.centroid_valid();
        Command.end('Calculated {0} spot centroids'.format(len(shoeboxes)))

        # Calculate the spot intensities
        Command.start('Calculating {0} spot intensities'.format(len(shoeboxes)))
        intensity = shoeboxes.summed_intensity_valid();
        Command.end('Calculated {0} spot intensities'.format(len(shoeboxes)))

        # Create the observations
        observed = flex.observation(shoeboxes.panels(), centroid, intensity)

        # Filter the reflections and select only the desired spots
        flags = self.filter_spots(observations=observed, shoeboxes=shoeboxes)
        observed = observed.select(flags)
        shoeboxes = shoeboxes.select(flags)

        # Return as a reflection list
        return ReflectionList(observed, shoeboxes)


class Filter(object):

    def __call__(self, flags, **kwargs):
        return self.run(self.check_flags(flags, **kwargs), **kwargs)

    def check_flags(self, flags, predictions=None, observations=None,
                    shoeboxes=None, **kwargs):
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


class NullFilter(Filter):

    def run(self, flags, **kwargs):
        return flags


class MinPixelsFilter(Filter):

    def __init__(self, num, code):
        self.code = code
        self.num = num

    def run(self, flags, observations=None, shoeboxes=None, **kwargs):

        # Get the number of mask values matching the code
        count = shoeboxes.count_mask_values(self.code)

        # Return the flags of those > the given number
        return flags.__and__(count > self.num)


class PeakCentroidDistanceFilter(Filter):

    def __init__(self, maxd):
        self.maxd = maxd

    def run(self, flags, observations=None, shoeboxes=None, **kwargs):

        # Get the peak locations and the centroids
        peak = shoeboxes.individual_peak_indices()
        cent = observations.centroids().px_positions()

        # Return the flags of those closer than the min distance
        return flags.__and__((peak - cent).norms() < self.max_d)


class CentroidResolutionFilter(Filter):

    def __init__(self, d_min, d_max):
        self.d_min = d_min
        self.d_max = d_max

    def run(self, flags, sweep=None, observations=None, **kwargs):

        # Get all the observation resolutions
        d = observations.resolutions(sweep.get_beam(), sweep.get_detector())

        # Return the flags of those in range
        return flags.__and__(d > self.d_min).__and__(d < self.d_max)


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
        # Configure the algorithm and wrap it up
        find_spots = SpotFinderFactory.configure_algorithm(params)
        return SpotFinder(
            find_spots=find_spots,
            filter_spots=NullFilter(),
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
        from dials.algorithms.peak_finding.spot_finder import FindSpots

        # Read in the lookup files
        gain_map = SpotFinderFactory.load_image(params.lookup.gain_map)
        dark_map = SpotFinderFactory.load_image(params.lookup.dark_map)

        # Create the threshold strategy
        threshold = SpotFinderFactory.configure_threshold(params, gain_map)

        # Setup the spot finder
        return FindSpots(threshold_image=threshold)

    @staticmethod
    def configure_threshold(params, gain_map):
        '''Get the threshold strategy'''
        from dials.algorithms.peak_finding.threshold \
            import UnimodalThresholdStrategy, XDSThresholdStrategy

        # Chose the strategy
        if params.spotfinder.threshold.algorithm == 'xds':
            return XDSThresholdStrategy(
                kernel_size=params.spotfinder.threshold.kernel_size,
                gain=gain_map,
                n_sigma_b=params.spotfinder.threshold.sigma_background,
                n_sigma_s=params.spotfinder.threshold.sigma_strong)

        elif params.spotfinder.threshold.algorithm == 'unimodal':
            return UnimodalThresholdStrategy()

        else:
            raise RuntimeError('Unknown threshold strategy')

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


class SpotFinderWrapperOld(object):
    ''' Temporary class while I tidy up spot finding '''

    def __init__(self, find_spots, scan_range):
        ''' Initialise '''
        self.find_spots = find_spots
        self.scan_range = scan_range

    def __call__(self, sweep):
        ''' Do the spot finding '''
        from dials.model.data import ReflectionList

        # Get list of scan ranges
        if not self.scan_range:
            scan_range = [(0, len(sweep))]
        else:
            scan_range = self.scan_range

        # Get spots from bits of scan
        observed = []
        for scan in scan_range:
            j0, j1 = scan
            assert(j1 > j0 and j0 >= 0 and j1 <= len(sweep))
            print '\nFinding spots in image {0} to {1}...'.format(j0, j1)

            # Calling spot finder algorithms
            observed.extend(self.find_spots(sweep[j0:j1]))

        # Return the spots in a reflection list
        return ReflectionList(observed)


class SpotFinderFactoryOld(object):
    ''' Factory class to create spot finders '''

    @staticmethod
    def from_parameters(params):
        ''' Given a set of parameters, construct the spot finder

        Params:
            params The input parameters

        Returns:
            The spot finder instance

        '''
        # Configure the algorithm and wrap it up
        find_spots = SpotFinderFactoryOld.configure_algorithm(params)
        return SpotFinderWrapperOld(find_spots, params.spotfinder.scan_range)

    @staticmethod
    def configure_algorithm(params):
        ''' Given a set of parameters, construct the spot finder

        Params:
            params The input parameters

        Returns:
            The spot finder instance

        '''
        if params.spotfinder.threshold.algorithm == 'lui':
            return SpotFinderFactoryOld.configure_lui_spotfinder(params)
        else:
            return SpotFinderFactoryOld.configure_general_spotfinder(params)

    @staticmethod
    def configure_lui_spotfinder(params):
        '''Configure the lui spotfinder algorithm.'''
        from dials.algorithms.peak_finding.spot_finder_lui import SpotFinderLui
        find_spots = SpotFinderLui(
            times = params.spotfinder.threshold.times,
            shift = params.spotfinder.threshold.shift,
            n_blocks_x = params.spotfinder.threshold.block_size[0],
            n_blocks_y = params.spotfinder.threshold.block_size[1],
            dimensions = params.spotfinder.threshold.dimensions)
        return find_spots

    @staticmethod
    def configure_general_spotfinder(params):
        '''Configure the general spot finder algorithms.'''
        from dials.util.command_line import Command
        from dials.algorithms.peak_finding.spot_finder import SpotFinderOld

        # Read in the lookup files
        gain_map = SpotFinderFactoryOld.load_image(params.lookup.gain_map)
        dark_map = SpotFinderFactoryOld.load_image(params.lookup.dark_map)

        # Create the threshold strategy
        threshold = SpotFinderFactoryOld.configure_threshold(params, gain_map)

        # Setup the spot finder
        return SpotFinderOld(
            min_spot_size=params.spotfinder.filter.min_spot_size,
            max_separation=params.spotfinder.filter.max_separation,
            d_min=params.spotfinder.filter.d_min,
            d_max=params.spotfinder.filter.d_max,
            threshold_strategy=threshold)

    @staticmethod
    def configure_threshold(params, gain_map):
        '''Get the threshold strategy'''
        from dials.algorithms.peak_finding.threshold \
            import UnimodalThresholdStrategy, XDSThresholdStrategy

        # Chose the strategy
        if params.spotfinder.threshold.algorithm == 'xds':
            return XDSThresholdStrategy(
                kernel_size=params.spotfinder.threshold.kernel_size,
                gain=gain_map,
                n_sigma_b=params.spotfinder.threshold.sigma_background,
                n_sigma_s=params.spotfinder.threshold.sigma_strong)

        elif params.spotfinder.threshold.algorithm == 'unimodal':
            return UnimodalThresholdStrategy()

        else:
            raise RuntimeError('Unknown threshold strategy')

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
