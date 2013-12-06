#!/usr/bin/env python
#
#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

"""This script has the single purpose of aligning a detector from a multipanel
sweep with a single panel detector, for a refinement test.

The multipanel detector has been created by the FormatCBFMiniPilatusSplit6M
class. It has 60 co-planar panels in the plane formed by its fast, slow axes.
The panel whose origin in that plane is the minimum coordinate in the basis
fast, slow is the 60th panel. Thus, to align with the single panel detector, the
origin of this panel must coincide with the origin of the single panel detector.
All panels must also share the same fast, slow directions as the single panel
detector.

The origin and fast, slow directions of the single panel detector are hard-coded
here. See dials_regression/refinement_test_data/i04_weak_data for the multi-
panel refinement test."""

from __future__ import division
import sys

if __name__ == '__main__':

  from optparse import OptionParser

  # Specify the command line options
  usage  = "usage: %prog [options] sweep.json"

  # Create an option parser
  parser = OptionParser(usage)

  parser.add_option('-o', '--output-file',
                  dest='output_file', type="string", default=None,
                  help='Destination filename for sweep with aligned detector')

  # Parse the arguments
  options, args = parser.parse_args()

  # Print help if not enough arguments specified
  if len(args) < 1:
    parser.print_help()
    sys.exit()

  from dials.model.serialize import load, dump

  # Try to load the models
  print 'Loading detector from {0}'.format(args[0])

  with open(args[0], 'r') as f: sweep = load.sweep(f)
  detector = sweep.get_detector()

  # calculate offsets
  from scitbx import matrix

  fast = matrix.col(detector[0].get_fast_axis()).normalize()
  slow = matrix.col(detector[0].get_slow_axis()).normalize()
  panel60_origin = matrix.col(detector[59].get_origin())

  fast_offsets = [(matrix.col(p.get_origin()) - panel60_origin).dot(fast) \
                  for p in detector]
  slow_offsets = [(matrix.col(p.get_origin()) - panel60_origin).dot(slow) \
                  for p in detector]

  ref_origin = matrix.col((-210.68270819092896,
                            205.69567063383158,
                           -263.7163617113081))
  ref_fast = matrix.col((0.9999974877443938,
                        -0.0015832096467848918,
                        -0.0015868056325671835))
  ref_slow = matrix.col((-0.0015900007754938183,
                         -0.9999895438652769,
                         -0.004287663425368496))

  for f, s, panel in zip(fast_offsets, slow_offsets, detector):

    new_origin = ref_origin + f * ref_fast + s * ref_slow
    panel.set_frame(ref_fast, ref_slow, new_origin)

  # Save the refined geometry to file
  print "detector aligned"
  if options.output_file:
    print 'Saving aligned detector geometry to {0}'.format(
        options.output_file)
    dump.sweep(sweep, open(options.output_file, 'w'))
