#!/usr/bin/env python
#
# detector_max_resolution.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

def print_detector_max_resolution(input_filename):
    """Read the given data file and print the resolution of the detector"""
    import dxtbx

    # Read the models from the input file
    print "Reading: \"{0}\"".format(input_filename)
    models = dxtbx.load(input_filename)
    beam = models.get_beam()
    detector = models.get_detector()
    assert(len(detector) == 1)

    # Print the resolution data
    print ''
    print detector
    print beam
    print ''
    print 'Beam centre: {0}'.format(detector[0].get_beam_centre(beam.get_s0()))
    print 'Max resolution at detector corners: {0}'.format(
        detector[0].get_max_resolution_at_corners(beam.get_s0()))
    print 'Max resolution for fully recorded elipse: {0}'.format(
        detector[0].get_max_resolution_ellipse(beam.get_s0()))

if __name__ == '__main__':
    import sys
    print_detector_max_resolution(sys.argv[1])
