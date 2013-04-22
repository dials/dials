#!/usr/bin/env python
#
# background_lookup.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.


class ScriptRunner(object):
    '''Class to run background lookup table script.'''

    def __init__(self, sweep_filenames):
        '''Setup the script.'''

        # Filename data
        self.sweep_filenames = sweep_filenames

    def __call__(self):
        '''Run the script.'''
        from dxtbx.sweep import SweepFactory
        from dials.util.command_line import Command, ProgressBar
        from dials.algorithms.background_lookup import ComputeDetectorNoise
        from dials.algorithms.background_lookup import ComputeBackgroundAndGain
        from scitbx.array_family import flex

        # Set the print output
        print 'Running {0}\n'.format(__file__)

        # Load the sweep
        Command.start('Loading sweep')
        sweep = SweepFactory.sweep(self.sweep_filenames)
        Command.end('Loaded sweep of {0} images.'.format(len(sweep)))

        # Any image pixels < 0 are set to 0 in the mask
        mask = sweep[0] >= 0

        # Create the objects for doing the calculations
        compute_noise = ComputeDetectorNoise()
        compute_lookup = ComputeBackgroundAndGain(mask)

        # Add images from the sweep to the calculations
        progress = ProgressBar('Calculating lookup')
        for i, image in enumerate(sweep):
            compute_noise.add(image)
            compute_lookup.add(image)
            progress.update(100.0 * (i + 1) / len(sweep))

        # Calculate the detector noise/gain/background
        noise = compute_noise.compute()
        gain = compute_lookup.gain()
        background = compute_lookup.background()
        progress.finished('Calculated lookup from {0} frames'.format(
            len(sweep)))

        # Display the gain
        from matplotlib import pylab, cm
        pylab.imshow(gain.as_numpy_array(), vmin=0, vmax=2, cmap=cm.Greys_r)
        ax = pylab.gca()
        ax.format_coord = lambda x, y: "x={0:.2f} y={1:.2f} value={2:.4f}"\
            .format(x, y, gain.as_numpy_array()[y, x])
        pylab.show()

        # Display the background
        pylab.imshow(background.as_numpy_array(), vmin=0, vmax=20,
            cmap=cm.Greys_r)
        ax = pylab.gca()
        ax.format_coord = lambda x, y: "x={0:.2f} y={1:.2f} value={2:.4f}"\
            .format(x, y, background.as_numpy_array()[y, x])
        pylab.show()


if __name__ == '__main__':

    from optparse import OptionParser

    # Specify the command line options
    usage  = "usage: %prog [options] /path/to/image/files "

    # Create an option parser
    parser = OptionParser(usage)

    # Parse the arguments
    options, args = parser.parse_args()

    # Print help if no arguments specified, otherwise call function
    if len(args) < 1:
        parser.print_help()
    else:
        # Initialise the script runner
        runner = ScriptRunner(sweep_filenames=args[0:])

        # Run the script
        runner()
