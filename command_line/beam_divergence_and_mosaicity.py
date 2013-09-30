from __future__ import division
#!/usr/bin/env python
#
# beam_divergence_and_mosaicity.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

class ScriptRunner(object):
    '''Class to run the beam divergence and mosaicity script.'''

    def __init__(self, spot_filename, sweep_filename):
        '''Setup the script.'''

        # Filename data
        self.spot_filename = spot_filename
        self.sweep_filename = sweep_filename

    def __call__(self):
        '''Run the script.'''
        from dxtbx.imageset import ImageSetFactory
        from dials.util.command_line import Command
        from dials.algorithms.integration.divergence_and_mosaicity \
            import BeamDivergenceAndMosaicity
        from dials.util.nexus import NexusFile
        import cPickle as pickle

        # Set the print output
        Command.indent = 4
        print 'Running {0}\n'.format(__file__)

        # Load the sweep
        Command.start('Loading sweep')
        sweep = ImageSetFactory.new(self.sweep_filename)
        assert(len(sweep) == 1)
        sweep = sweep[0]
        Command.end('Loaded sweep of {0} images.'.format(len(sweep)))

        # Heading for file read bit
        print 'Reading datafiles...'
        Command.start('Loading reflections')
        reflections = pickle.load(open(self.spot_filename, 'rb'))
        Command.end('Loaded {0} reflections'.format(len(reflections)))
        # Get the reflection list
#        Command.start('Reading reflection file.')
#        nexus_handle = NexusFile(self.spot_filename, 'r')
#        predicted = nexus_handle.get_reflections()
#        nexus_handle.close()
#        Command.end('Read {0} predicted reflections.'.format(len(predicted)))

        # Setup the beam divergence and mosaicity calculator
        compute_sigma = BeamDivergenceAndMosaicity(sweep)

        # Calculate the e.s.d of beam divergence and mosaicity
        compute_sigma(reflections)


if __name__ == '__main__':

    from optparse import OptionParser

    # Specify the command line options
    usage  = "usage: %prog [options] " \
             "/path/to/reflection_file.pkl " \
             "/path/to/image_files"

    # Create an option parser
    parser = OptionParser(usage)

    # Parse the arguments
    options, args = parser.parse_args()

    # Print help if no arguments specified, otherwise call function
    if len(args) < 2:
        parser.print_help()
    else:
        # Initialise the script runner
        runner = ScriptRunner(spot_filename=args[0], sweep_filename=args[1:])

        # Run the script
        runner()
