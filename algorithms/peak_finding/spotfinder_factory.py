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

            # Filter the spots
            spots = self.filter_spots(spots)

            # Add the spots to the list
            spots_all.extend(spots)

        # Calculate the centroids
        Command.start('Calculating {0} centroids'.format(len(spots_all)))
        cpos, cvar, cerr, ctot = self.centroid(spots_all)
        Command.end('Calculated {0} centroids'.format(len(spots_all)))

        # Return the spots in a reflection list
        Command.start('Creating reflection list')
        rlist = self.reflection_list(spots_all, cpos, cvar, cerr, ctot)
        Command.end('Created list of {0} reflections'.format(len(rlist)))
        return rlist

    def centroid(self, spots):
        '''Calculate the spot centroids.

        Params:
            spots The list of spots

        Returns:
            (centroid position, centroid variance)

        '''
        from dials.algorithms.image.centroid import CentroidMaskedImage3d
        from scitbx.array_family import flex

        # Initialise arrays
        centroid_pos = flex.vec3_double()
        centroid_var = flex.vec3_double()
        centroid_err = flex.vec3_double()
        centroid_tot = flex.double()

        # Loop through each spot
        for s in spots:

            # Find the spot centroid
            centroid = CentroidMaskedImage3d(s.data, s.mask)
            pos = centroid.mean()
            pos = pos[0] + s.bbox[0], pos[1] + s.bbox[2], pos[2] + s.bbox[4]
            centroid_pos.append(pos)
            centroid_tot.append(centroid.sum_pixels())
            try:
                centroid_var.append(centroid.unbiased_variance())
                centroid_err.append(centroid.unbiased_standard_error_sq())
            except RuntimeError:
                centroid_var.append((0.5, 0.5, 0.5))
                centroid_err.append((0.5, 0.5, 0.5))

        # Return the centroid and variance
        return centroid_pos, centroid_var, centroid_err, centroid_tot

    def reflection_list(self, spots, cpos, cvar, cerr, ctot):
        '''Create a reflection list from the spot data.

        Params:
            spots The spot list
            cpos The centroid position
            cvar The centroid variance
            cerr The centroid error
            ctot The centroid total counts

        Returns:
            A list of reflections

        '''
        from dials.model.data import Reflection, ReflectionList
        from dials.algorithms import shoebox

        # Ensure the lengths are ok
        assert(len(spots) > 0)
        assert(len(spots) == len(cpos))
        assert(len(spots) == len(cvar))
        assert(len(spots) == len(cerr))
        assert(len(spots) == len(ctot))

        # Create the reflection list
        rlist = ReflectionList(len(spots))
        for i in range(len(spots)):

            # Set the shoebox info
            rlist[i].bounding_box = spots[i].bbox
            rlist[i].shoebox = spots[i].data
            rlist[i].shoebox_mask = spots[i].mask

            # Set the centroid and intensity info
            rlist[i].centroid_position = cpos[i]
            rlist[i].centroid_variance = cerr[i]
            rlist[i].centroid_sq_width = cvar[i]
            rlist[i].intensity = ctot[i]

        # Return the reflection list
        return rlist


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
