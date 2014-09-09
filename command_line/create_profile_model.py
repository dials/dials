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
  from dials.util.command_line import Importer
  from optparse import OptionParser
  from dials.algorithms.profile_model.profile_model import ProfileModelList
  from math import pi
  from dials.util.command_line import Command
  from dials.array_family import flex

  # The script usage
  usage  = "usage: %prog [options] experiments.json spots.pickle"
  parser = OptionParser(usage=usage)

  # Output filename option
  parser.add_option(
    '-o', '--output-filename',
    dest = 'output_filename',
    type = 'string', default = 'profile.phil',
    help = 'Set the filename for profile parameters')

  # Parse the command line
  options, args = parser.parse_args()

  # Import the items
  importer = Importer(args, check_format=False)
  experiments = importer.experiments
  if experiments is None or len(experiments) != 1:
    parser.print_help()
    exit(0)
  if importer.reflections is None or len(importer.reflections) != 1:
    parser.print_help()
    exit(0)
  reflections = importer.reflections[0]

  from dials.array_family import flex
  Command.start('Removing invalid coordinates')
  xyz = reflections['xyzcal.mm']
  mask = flex.bool([x == (0, 0, 0) for x in xyz])
  reflections.del_selected(mask)
  Command.end('Removed invalid coordinates, %d remaining' % len(reflections))

  # Create the profile model
  profile_model = ProfileModelList.compute(experiments, reflections)
  for model in profile_model:
    sigma_b = model.sigma_b(deg=True)
    sigma_m = model.sigma_m(deg=True)
    print 'Sigma B: %f' % sigma_b
    print 'Sigma M: %f' % sigma_m

  # Wrtie the parameters
  filename = options.output_filename
  Command.start("Writing profile model to %s" % filename)
  with open(filename, "w") as outfile:
    outfile.write(profile_model.dump())
  Command.end("Wrote profile model to %s" % filename)
