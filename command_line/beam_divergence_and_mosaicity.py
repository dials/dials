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

    def __init__(self, spot_filename, sweep_filenames, min_spot_size = 6,
        max_pc_separation = 2, max_bc_separation = 2):
        '''Setup the script.'''

        # Filename data
        self.spot_filename = spot_filename
        self.sweep_filenames = sweep_filenames

        # Set spot finder parameters
        self.min_spot_size = min_spot_size
        self.max_pc_separation = max_pc_separation

        # Set beam divergence/mosaicity parameters
        self.max_bc_separation = max_bc_separation

    def print_params(self):
        '''Print the parameters.'''
        print 'Parameters...'
        print '    Minimum spot size = {0}.'.format(self.min_spot_size)
        print '    Maximum peak-centroid distance = {0}'\
            .format(self.max_pc_separation)
        print '    Maximum bragg peak-centroid distance = {0}'\
            .format(self.max_bc_separation)
        print ''
        
    def __call__(self):
        '''Run the script.'''
        from dxtbx.sweep import SweepFactory
        from dials.util.command_line import Command
        from dials.algorithms.peak_finding.spot_finder import SpotFinder
        from dials.algorithms.integration.divergence_and_mosaicity \
            import BeamDivergenceAndMosaicity
        from dials.util.nexus import NexusFile
        
        # Set the print output
        Command.indent = 4
        print 'Running {0}\n'.format(__file__)
        self.print_params()
        
        # Heading for file read bit       
        print 'Reading datafiles...'
        
        # Load the sweep
        Command.start('Loading sweep')
        sweep = SweepFactory.sweep(self.sweep_filenames)
        Command.end('Loaded sweep of {0} images.'.format(len(sweep)))

        # Get the reflection list
        Command.start('Reading reflection file.')
        nexus_handle = NexusFile(self.spot_filename, 'r')
        predicted = nexus_handle.get_reflections()
        nexus_handle.close()
        Command.end('Read {0} predicted reflections.'.format(len(predicted)))

        # Setup the algorithms
        find_spots = SpotFinder(
            min_spot_size=self.min_spot_size,
            max_separation=self.max_pc_separation)
        compute_sigma = BeamDivergenceAndMosaicity(
            sweep,
            max_separation=self.max_bc_separation)

        # Find some spots in the sweep
        observed = find_spots(sweep)

        # Calculate the e.s.d of beam divergence and mosaicity
        compute_sigma(observed, predicted)


if __name__ == '__main__':

    from optparse import OptionParser

    # Specify the command line options
    usage  = "usage: %prog [options] /path/to/reflection_file.h5 "
    usage += "/path/to/image.cbf"

    # Create an option parser
    parser = OptionParser(usage)

    # Add some command line options
    parser.add_option(
        '-s', '--min-spot-size',
        dest = 'min_spot_size', type = "int", default = 6,
        help = 'Minimum pixels in spot')
    parser.add_option(
        '-p', '--max-pc-separation',
        dest = 'max_pc_separation', type = "int", default = 2,
        help = 'Maximum peak-centroid distance')
    parser.add_option(
        '-b', '--max-bc-separation',
        dest = 'max_bc_separation', type = "int", default = 2,
        help = 'Maximum bragg peak-centroid distance')

    # Parse the arguments
    options, args = parser.parse_args()

    # Print help if no arguments specified, otherwise call function
    if len(args) < 2:
        print parser.print_help()
    else:
        # Initialise the script runner
        runner = ScriptRunner(
            spot_filename=args[0],
            sweep_filenames=args[1:],
            min_spot_size=options.min_spot_size,
            max_pc_separation=options.max_pc_separation,
            max_bc_separation=options.max_bc_separation)

        # Run the script
        runner()

