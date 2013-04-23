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

class ScriptRunner(object):
    '''Class to run script.'''

    def __init__(self, sweep_filenames, **kwargs):
        '''Setup the script.'''

        # Filename data
        self.sweep_filenames = sweep_filenames

        # Set the options
        self.threshold         = kwargs['threshold']
        self.output_file       = kwargs['output_file']

        # Threshold options
        self.dark_current_file = kwargs['dark_current_file']
        self.gain_map_file     = kwargs['gain_map_file']
        self.sigma_background  = kwargs['sigma_background']
        self.sigma_strong      = kwargs['sigma_strong']
        self.kernel_size       = kwargs['kernel_size']

        # Filter options
        self.min_spot_size     = kwargs['min_spot_size']
        self.max_pc_separation = kwargs['max_pc_separation']

        # Set some other stuff
        self.gain_map = None
        self.dark_current = None
        self.mask = None

    def __call__(self):
        '''Run the script.'''
        from scitbx.array_family import flex
        from dxtbx.sweep import SweepFactory
        from dials.util.command_line import Command
        from dials.algorithms.peak_finding.spot_finder import SpotFinder
        from dials.algorithms.peak_finding.threshold import XDSThresholdStrategy

        # Set the print output
        Command.indent = 4
        print 'Running {0}\n'.format('spotfinder')

        # Heading for file read bit
        print 'Reading datafiles...'

        # Load the sweep
        Command.start('Loading sweep')
        sweep = SweepFactory.sweep(self.sweep_filenames)
        Command.end('Loaded sweep of {0} images.'.format(len(sweep)))

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
            self.gain = handle.get_data(dim)
            Command.end('Read gain map file')

        # Create the mask
        self.mask = sweep[0] >= 0

        # Setup the spot finder
        find_spots = self._spotfinder(sweep)

        # Find some spots in the sweep
        observed = find_spots(sweep)

        # Save the reflections
        if self.output_file:
            import pickle
            Command.start('Saving spots to {0}'.format(self.output_file))
            pickle.dump(observed, open(self.output_file, 'wb'))
            Command.end('Saved spots to {0}'.format(self.output_file))

    def _spotfinder(self, sweep):
        '''Get the spot finder'''
        from dials.algorithms.peak_finding.spot_finder import SpotFinder
        return SpotFinder(
            min_spot_size = self.min_spot_size,
            max_separation = self.max_pc_separation,
            threshold_strategy = self._threshold_strategy(sweep))

    def _threshold_strategy(self, sweep):
        '''Get the threshold strategy'''
        from dials.algorithms.peak_finding.threshold \
            import UnimodalThresholdStrategy, XDSThresholdStrategy

        # Chose the strategy
        if self.threshold == 'xds':
            return XDSThresholdStrategy(
                kernel_size = self.kernel_size,
                gain = self.gain_map,
                mask = self.mask,
                n_sigma_b = self.sigma_background,
                n_sigma_s = self.sigma_strong)
        elif self.threshold == 'unimodal':
            trusted_range = sweep.get_detector().get_trusted_range()
            return UnimodalThresholdStrategy(trusted_range=trusted_range)
        else:
            raise RuntimeError('Unknown threshold strategy')


if __name__ == '__main__':

    from optparse import OptionParser, OptionGroup, IndentedHelpFormatter

    # Specify the command line options
    usage  = "usage: %prog [options] /path/to/image/files"

    # Create an option parser
    parser = OptionParser(usage)

    # Add algorithm options
    parser.add_option(
        '-t', '--threshold',
        dest = 'threshold',
        type = 'choice', choices = ['xds', 'unimodal'], default = 'xds',
        help = 'The threshold algorithm to use (default = %default)')
    parser.add_option(
        '-o', '--output-file',
        dest = 'output_file',
        type = 'string', default = None,
        help = 'Select a file to save the spots.')

    # Create a group for threshold options
    threshold_group = OptionGroup(
        parser, 'Threshold Options',
        'Options affecting threshold algorithms')
    threshold_group.add_option(
        '--dark-current-file',
        dest = 'dark_current_file', type = "string", default = None,
        help = 'Dark current map filename')
    threshold_group.add_option(
        '--gain-map-file',
        dest = 'gain_map_file', type = "string", default = None,
        help = 'Gain map filename')
    threshold_group.add_option(
        '--sigma-background',
        dest = 'sigma_background', type = 'float', default = 6.0,
        help = '(var/mean) > gain + n_sigma*gain*sqrt(2/(n - 1))')
    threshold_group.add_option(
        '--sigma-strong',
        dest = 'sigma_strong', type = 'float', default = 3.0,
        help = 'pixel > mean + n_sigma*sdev (used by: xds) (default: %default)')
    threshold_group.add_option(
        '--kernel-size',
        dest = 'kernel_size', type = 'int', nargs=2, default = (3, 3),
        help = 'Local window size (2 * s + 1) centred on pixel')

    # Create group for filter options
    filter_group = OptionGroup(
        parser, 'Filter options',
        'Options affecting filter algorithms')
    filter_group.add_option(
        '--min-spot-size',
        dest = 'min_spot_size', type = "int", default = 6,
        help = 'Minimum pixels in spot')
    filter_group.add_option(
        '--max-pc-separation',
        dest = 'max_pc_separation', type = "int", default = 2,
        help = 'Maximum peak-centroid distance')

    # Add the related groups of options
    parser.add_option_group(threshold_group)
    parser.add_option_group(filter_group)

    # Parse the arguments
    options, args = parser.parse_args()

    # Print help if no arguments specified, otherwise call function
    if len(args) < 1:
        parser.print_help()
    else:

        # Initialise the script runner
        runner = ScriptRunner(sweep_filenames=args, **options.__dict__)

        # Run the script
        runner()
