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

  #xy = reflections.compute_ray_intersections(experiments[0].detector)
  #from dials.algorithms.spot_prediction import RotationAngles

  #s0 = experiments[0].beam.get_s0()
  #m2 = experiments[0].goniometer.get_rotation_axis()
  #compute_angles = RotationAngles(s0, m2)
  #phis = [compute_angles(h, ub) for h in reflections['miller_index']]
  #obsphis = reflections['
  reflections['xyzcal.mm'] = reflections['xyzobs.mm.value']

  # Create the profile model
  profile_model = ProfileModel(experiments[0], reflections)
  print 'Sigma B: %f' % profile_model.sigma_b()
  print 'Sigma M: %f' % profile_model.sigma_m()
