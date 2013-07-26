#!/usr/bin/env python
#
# spotfinder.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME dials.spotfinder
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

class ScriptRunner(object):
    '''Class to run script.'''

    def __init__(self, sweep_filenames, output_file, params):
        '''Setup the script.'''

        # Filename data
        self.sweep_filenames = sweep_filenames
        self.output_file = output_file

        # Set the options
        self.algorithm = params.threshold.algorithm
        self.scan_range = params.scan_range

        # Threshold options
        self.dark_current_file = params.threshold.noise_file
        self.gain_map_file = params.threshold.gain_file
        self.sigma_background = params.threshold.sigma_background
        self.sigma_strong = params.threshold.sigma_strong
        self.kernel_size = params.threshold.kernel_size
        self.times = params.threshold.times
        self.shift = params.threshold.shift
        self.block_size = params.threshold.block_size
        self.dimensions = params.threshold.dimensions

        # Filter options
        self.min_spot_size = params.filter.min_spot_size
        self.max_pc_separation = params.filter.max_separation
        self.d_min = params.filter.d_min
        self.d_max = params.filter.d_max

        # Set some other stuff
        self.gain_map = None
        self.dark_current = None
        self.mask = None

        self.params = params

    def __call__(self):
        '''Run the script.'''
        from dxtbx.imageset import ImageSetFactory
        from dials.util.command_line import Command
        from dials.algorithms.peak_finding.spot_finder import SpotFinder
        from dials.algorithms.peak_finding.threshold import XDSThresholdStrategy
        from dials.model.data import ReflectionList

        # Set the print output
        Command.indent = 4
        print 'Running {0}\n'.format('spotfinder')

        # Heading for file read bit
        print 'Reading datafiles...'

        if len(self.sweep_filenames) == 1:
            raise RuntimeError("spotfinding currently requires more than one image.")

        # Load the sweep
        Command.start('Loading sweep')
        sweep = ImageSetFactory.new(self.sweep_filenames)
        assert(len(sweep) == 1)
        sweep = sweep[0]
        Command.end('Loaded sweep of {0} images.'.format(len(sweep)))

        # Get list of scan ranges
        if len(self.scan_range) > 0:
            scan_range = self.scan_range
        else:
            scan_range = [(0, len(sweep))]

        # Get spots from bits of scan
        observed = []
        for scan in scan_range:
            j0, j1 = scan
            assert(j1 > j0 and j0 >= 0 and j1 <= len(sweep))
            print '\nFinding spots in image {0} to {1}...'.format(j0, j1)

            # Calling spot finder algorithms
            if self.algorithm == 'lui':
                observed.extend(self._spotfinder_lui(sweep[j0:j1]))
            else:
                observed.extend(self._spotfinder_xds(sweep[j0:j1]))

        # Store the reflections
        self.reflections = observed

        if self.params.image_viewer:
            self.view()

        # Save the reflections
        if self.output_file:
            import cPickle as pickle
            Command.start('Saving spots to {0}'.format(self.output_file))
            pickle.dump(ReflectionList(observed), open(self.output_file, 'wb'),
                pickle.HIGHEST_PROTOCOL)
            Command.end('Saved spots to {0}'.format(self.output_file))

    def _spotfinder_lui(self, sweep):
        '''Call the lui spotfinder algorithm.'''
        from dials.algorithms.peak_finding.spot_finder_lui import SpotFinderLui
        find_spots = SpotFinderLui(times = self.times, shift = self.shift,
          n_blocks_x = self.block_size[0], n_blocks_y = self.block_size[1],
          dimensions = self.dimensions)
        return find_spots(sweep)

    def _spotfinder_xds(self, sweep):
        '''Call the XDS spot finder algorithms.'''
        from dials.util.command_line import Command

        # Read the dark current file
        if self.dark_current_file:
            Command.start('Reading dark current file')
            from iotbx.xds import blank_cbf
            handle = blank_cbf.reader()
            handle.read_file(self.dark_current_file)
            self.dark_current = handle.get_data()
            Command.end('Read dark current file')

        # Read the gain map file
        if self.gain_map_file:
            from iotbx.xds import gain_cbf
            Command.start('Reading gain map file')
            handle = gain_cbf.reader()
            handle.read_file(self.gain_map_file)
            dim = sweep.get_detector().get_image_size()[::-1]
            self.gain_map = handle.get_data(dim)
            Command.end('Read gain map file')

        # Create the mask
        self.mask = sweep[0] >= 0

        # Setup the spot finder
        find_spots = self._spotfinder(sweep)

        # Find some spots in the sweep
        return find_spots(sweep)

    def _spotfinder(self, sweep):
        '''Get the spot finder'''
        from dials.algorithms.peak_finding.spot_finder import SpotFinder
        return SpotFinder(
            min_spot_size = self.min_spot_size,
            max_separation = self.max_pc_separation,
            d_min = self.d_min,
            d_max = self.d_max,
            threshold_strategy = self._threshold_strategy(sweep))

    def _threshold_strategy(self, sweep):
        '''Get the threshold strategy'''
        from dials.algorithms.peak_finding.threshold \
            import UnimodalThresholdStrategy, XDSThresholdStrategy

        # Chose the strategy
        if self.algorithm == 'xds':
            return XDSThresholdStrategy(
                kernel_size = self.kernel_size,
                gain = self.gain_map,
                mask = self.mask,
                n_sigma_b = self.sigma_background,
                n_sigma_s = self.sigma_strong)
        elif self.algorithm == 'unimodal':
            trusted_range = sweep.get_detector().get_trusted_range()
            return UnimodalThresholdStrategy(trusted_range = trusted_range)
        else:
            raise RuntimeError('Unknown threshold strategy')

    def view(self):
        from dials.util.spotfinder_wrap import spot_wrapper
        spot_wrapper(working_phil=self.params).display(
            sweep_filenames=self.sweep_filenames, reflections=self.reflections)

if __name__ == '__main__':

    from dials.util.options import OptionParser

    # Specify the command line options
    usage = "usage: %prog [options] /path/to/image/files"

    # Create an option parser
    parser = OptionParser(home_scope='spotfinder', usage=usage)

    # Add options
    parser.add_option(
        '-o', '--output-file',
        dest = 'output_file',
        type = 'string', default = None,
        help = 'Select a file to save the spots.')

    # Parse the arguments
    params, options, args = parser.parse_args()

    # Print help if no arguments specified, otherwise call function
    if len(args) < 1:
        parser.print_help()
    else:

        # Initialise the script runner
        runner = ScriptRunner(
            sweep_filenames = args,
            output_file=options.output_file,
            params=params.spotfinder)

        # Run the script
        runner()
