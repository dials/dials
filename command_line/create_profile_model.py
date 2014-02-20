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
  from math import pi
  from dials.util.command_line import Command

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
  Command.start('Importing Data')
  importer = Importer(args)
  experiments = importer.experiments
  assert(len(experiments) == 1)
  assert(len(importer.reflections) == 1)
  reflections = importer.reflections[0]
  Command.end('Imported %d reflections' % len(reflections))

  from dials.array_family import flex
  Command.start('Removing invalid coordinates')
  xyz = reflections['xyzcal.mm']
  mask = flex.bool([x == (0, 0, 0) for x in xyz])
  reflections.del_selected(mask)
  Command.end('Removed invalid coordinates, %d remaining' % len(reflections))

  # Create the profile model
  profile_model = ProfileModel(experiments[0], reflections)
  sigma_b = profile_model.sigma_b() * 180.0 / pi
  sigma_m = profile_model.sigma_m() * 180.0 / pi
  print 'Sigma B: %f' % sigma_b
  print 'Sigma M: %f' % sigma_m

  # Write the parameters
  from dials.framework.registry import Registry
  registry = Registry()

  # Get the parameters
  params = registry.config().params()
  params.shoebox.sigma_b = sigma_b
  params.shoebox.sigma_m = sigma_m

  # Get the diff phil to save
  master_phil = registry.config().phil()
  modified_phil = master_phil.format(python_object=params)
  diff_phil = master_phil.fetch_diff(source=modified_phil)

  # Wrtie the parameters
  filename = options.output_filename
  Command.start("Writing profile model to %s" % filename)
  with open(filename, "w") as outfile:
    outfile.write(diff_phil.as_str())
  Command.end("Wrote profile model to %s" % filename)
