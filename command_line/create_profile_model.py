#!/usr/bin/env python
#
# dials.create_profile_model.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division


if __name__ == '__main__':
  from optparse import OptionParser
  from dials.util.command_line import Importer
  from dials.algorithms.profile_model.profile_model import ProfileModel

  # The script usage
  usage  = "usage: %prog [options] experiments.json spots.pickle"
  parser = OptionParser(usage=usage)

  # Parse the command line
  options, args = parser.parse_args()

  # Import the items
  importer = Importer(args)

  # Get the experiments
  experiments = importer.experiments
  assert(len(experiments) == 1)

  # Get the reflections
  assert(len(importer.reflections) == 1)
  reflections = importer.reflections[0]
  print 'Imported %d reflections' % len(reflections)

  # Create the profile model
  profile_model = ProfileModel(experiments[0], reflections)
  print 'Sigma B: %f' % profile_model.sigma_b()
  print 'Sigma M: %f' % profile_model.sigma_m()
