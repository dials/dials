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

  from dials.array_family import flex
  xyz = reflections['xyzcal.mm']
  mask = flex.bool([x == (0, 0, 0) for x in xyz])
  reflections.del_selected(mask)
  reflections_all = reflections
  print 'Have %d reflections remaining' % len(reflections)

  sigma_b = []
  sigma_m = []
  for i in range(1):
    import random
    #index = random.sample(range(len(reflections_all)), 5000)
    #reflections = reflections_all.select(flex.size_t(index))

    # Create the profile model
    profile_model = ProfileModel(experiments[0], reflections)
    print 'Sigma B: %f' % profile_model.sigma_b()
    print 'Sigma M: %f' % profile_model.sigma_m()
    sigma_b.append(profile_model.sigma_b())
    sigma_m.append(profile_model.sigma_m())

  #from matplotlib import pylab
  #pylab.subplot(121)
  #pylab.hist(sigma_b)
  #pylab.subplot(122)
  #pylab.hist(sigma_m)
  #pylab.show()
