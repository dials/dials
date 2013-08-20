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

        # Set the spot finding function
        assert(find_spots != None)
        self.find_spots = find_spots

        # Set the filter function
        if filter_spots == None:
            self.filter_spots = lambda x: x
        else:
            self.filter_spots = filter_spots

        # Set the scan range
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
        spots_all = []
        for scan in scan_range:
            j0, j1 = scan
            assert(j1 > j0 and j0 >= 0 and j1 <= len(sweep))
            print '\nFinding spots in image {0} to {1}...'.format(j0, j1)

            # Find the spots
            spots = self.find_spots(sweep[j0:j1])

            # Filter the spots
            spots = self.filter_spots(spots)

            # Add the spots to the list
            spots_all.extend(spots)

        # Return the spots in a reflection list
        return spots_all


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
        return SpotFinder(find_spots=find_spots,
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
