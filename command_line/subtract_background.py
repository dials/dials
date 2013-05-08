#
# dials.command_line.subtract_background.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

class ScriptRunner(object):
    '''Script to subtract background from reflections.'''

    def __init__(self, params, options, args):
        '''Initialise the script.'''

        # Get the input and output filenames
        self.input_file = args[0]
        self.output_file = options.output_file

        # Get the parameters
        self.algorithm = params.background.algorithm
        self.min_pixels = params.background.min_pixels
        self.num_sigma = params.background.sigma_background

    def __call__(self):
        '''Subtract the background from the reflections.'''

        from dials.util.command_line import Command
        from dials.model.data import Reflection, ReflectionList
        import pickle

        # Read in the reflections from the inputfile
        Command.start('Reading reflections from {0}'.format(self.input_file))
        reflections = pickle.load(open(self.input_file, 'rb'))
        Command.end('Read {0} reflections from {1}'.format(
            len(reflections), self.input_file))

        # Create the subtraction algorithm
        subtract_background = self.get_algorithm()

        # Subtract the background from each reflection
        Command.start('Subtracting background from {0} reflections'.format(
            len(reflections)))
        valid = subtract_background(reflections)
        Command.end('Subtracted background from {0} reflections'.format(
            len(valid)))

        # Save the reflections
        if self.output_file:
            Command.start('Saving reflections to {0}'.format(self.output_file))
            pickle.dump(reflections, open(self.output_file, 'wb'))
            Command.end('Saved reflections to {0}'.format(self.output_file))

    def get_algorithm(self):
        '''Get the background subtraction algorithm.'''
        from dials.algorithms.background.subtraction import BackgroundSubtractor

        from dials.algorithms.background.flat_background_subtraction import FlatBackgroundSubtraction
        # Choose the algorithm object
        if self.algorithm == 'xds':
            algorithm = BackgroundSubtractor(
                min_pixels = self.min_pixels,
                num_sigma = self.num_sigma)
        elif self.algorithm == 'flat':
            algorithm = FlatBackgroundSubtraction()
        elif self.algorithm == 'curved':
            pass

        # Return the algorithm object
        return algorithm

if __name__ == '__main__':

    from dials.util.options import OptionParser

    # Specify the command line options
    usage = "usage: %prog [options] /path/to/reflection/file"

    # Create an option parser
    parser = OptionParser(home_scope = 'background', usage = usage)

    # Add options
    parser.add_option(
        '-o', '--output-file',
        dest = 'output_file',
        type = 'string', default = None,
        help = 'Select an outout file.')

    # Parse the arguments
    params, options, args = parser.parse_args()

    # Print help if no arguments specified, otherwise call function
    if len(args) < 1:
        parser.print_help()
    else:

        # Initialise the script runner
        run = ScriptRunner(params, options, args)
        run()
