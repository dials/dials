#!/usr/bin/env python
#
# print_reflection_stats.py
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

    def __init__(self, filename):
        '''Setup the script.'''

        # Filename data
        self.filename = filename

    def __call__(self):
        '''Run the script.'''
        import pickle
        from scitbx.array_family import flex
        from dials.model.data import ReflectionList
        from math import sqrt

        # Load the reflection list
        rlist = pickle.load(open(self.filename, 'r'))

        # Create the format string
        fmt  = '{0:>8,.2f} {1:>8,.2f} {2:>8,.2f}'
        fmt += '{3:>8,.2f} {4:>8,.2f} {5:>8,.2f}'
        fmt += '{6:>8}'

        # Loop through all the reflections
        for r in rlist:

            # Get data for reflection
            x0, x1, y0, y1, z0, z1 = r.bounding_box
            xc, yc, zc = r.centroid_position
            vx, vy, vz = r.centroid_variance
            shoebox = r.shoebox

            # Get distance between centroid and box centre
            dx = xc - ((x0 + x1) / 2.0)
            dy = yc - ((y0 + y1) / 2.0)
            dz = zc - ((z0 + z1) / 2.0)

            # Total intensity
            itot = flex.sum(shoebox)

            # Calculate ratio of difference to sdev
            dx_by_sdev = dx / sqrt(vx)
            dy_by_sdev = dy / sqrt(vy)
            dz_by_sdev = dz / sqrt(vz)

            # Print to stdout. Try is to allow paging with less
            try:
                print fmt.format(dx, dy, dz,
                                 dx_by_sdev, dy_by_sdev, dz_by_sdev,
                                 itot)
            except IOError:
                break

if __name__ == '__main__':

    from optparse import OptionParser

    # Specify the command line options
    usage  = "usage: %prog [options] /path/to/reflections.pkl"

    # Create an option parser
    parser = OptionParser(usage)

    # Parse the arguments
    options, args = parser.parse_args()

    # Print help if no arguments specified, otherwise call function
    if len(args) < 1:
        parser.print_help()
    else:
        # Initialise the script runner
        runner = ScriptRunner(filename=args[0])

        # Run the script
        runner()
