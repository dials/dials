#!/usr/bin/env python
#
# match_observed_w_predicted.py
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

    def __init__(self, observed_file, predicted_file, **kwargs):
        '''Setup the script.'''

        # Filename data
        self.observed_file  = observed_file
        self.predicted_file = predicted_file
        self.output_file    = kwargs['output_file']

        # Set the options
        self.max_separation = kwargs['max_separation']

    def __call__(self):
        '''Run the script.'''
        from dials.model.data import ReflectionList
        from dials.algorithms.peak_finding.spot_matcher import SpotMatcher
        from dials.util.command_line import Command
        import pickle

        # Set the print output
        Command.indent = 4
        print 'Running {0}\n'.format('match_observed_w_predicted')

        # Heading for file read bit
        print 'Reading datafiles...'
        Command.start('Loading observed reflections')
        observed = pickle.load(open(self.observed_file, 'rb'))
        Command.end('Loaded {0} observed reflections'.format(len(observed)))

        Command.start('Loading predicted reflections')
        predicted = pickle.load(open(self.predicted_file, 'rb'))
        Command.end('Loaded {0} predicted reflections'.format(len(predicted)))

        print 'Processing reflections...'
        Command.start('Matching observed with predicted spots')
        match_spots = SpotMatcher(max_separation=self.max_separation)
        matched = match_spots(observed, predicted)
        Command.end('Matched {0} observed with predicted spots'
            .format(len(matched)))

        # Save the reflections
        if self.output_file:
            print 'Writing datafiles...'
            Command.start('Saving spots to {0}'.format(self.output_file))
            pickle.dump(matched, open(self.output_file, 'wb'))
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
    usage  = "usage: %prog [options] "
    usage += "/path/to/observed.pkl "
    usage += "/path/to/predicted.pkl"

    # Create an option parser
    parser = OptionParser(usage)

    # Add algorithm options
    parser.add_option(
        '-o', '--output-file',
        dest = 'output_file',
        type = 'string', default = None,
        help = 'Select a file to save the matched spots.')

    parser.add_option(
        '--max-separation',
        dest = 'max_separation', type = "int", default = 2,
        help = 'Maximum bragg peak-centroid distance')

    # Parse the arguments
    options, args = parser.parse_args()

    # Print help if no arguments specified, otherwise call function
    if len(args) < 2:
        parser.print_help()
    else:

        # Initialise the script runner
        runner = ScriptRunner(
            observed_file=args[0],
            predicted_file=args[1],
            **options.__dict__)

        # Run the script
        runner()
