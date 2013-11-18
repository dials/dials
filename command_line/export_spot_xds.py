#!/usr/bin/env python
#
# export_spot_xds.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: Richard Gildea
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

class ScriptRunner(object):
  '''Class to run script.'''

  def __init__(self, reflections_filename, output_filename):
    '''Setup the script.'''

    # Filename data
    self.reflections_filename = reflections_filename
    self.output_filename = output_filename

  def __call__(self):
    '''Run the script.'''
    import cPickle as pickle
    from iotbx.xds import spot_xds
    from dials.model.data import ReflectionList # import dependency
    from dials.util.command_line import Command

    # Read the pickle file
    Command.start('Reading reflection file.')
    with open(self.reflections_filename, 'rb') as f:
      reflections = pickle.load(f)

    Command.end('Read {0} spots from reflection file.'.format(len(reflections)))

    # Save the reflection list
    if self.output_filename != None:
      Command.start('Saving reflections to {0}'.format(
          self.output_filename))

      centroids = []
      intensities = []
      miller_indices = []
      miller_indices_excluding_zero = []

      for refl in reflections:
        centroids.append(refl.centroid_position)
        intensities.append(refl.intensity)
        miller_indices_excluding_zero.append(refl.miller_index)
        if refl.miller_index != (0,0,0):
          miller_indices.append(refl.miller_index)

      if len(miller_indices) == 0:
        miller_indices = None
      xds_writer = spot_xds.writer(centroids=centroids,
                                   intensities=intensities,
                                   miller_indices=miller_indices)
      xds_writer.write_file(filename=self.output_filename)

      Command.end('Saved reflections to {0}'.format(
          self.output_filename))

if __name__ == '__main__':

  from optparse import OptionParser

  # Specify the command line options
  usage  = "usage: %prog [options] " \
           "/path/to/reflections.pickle "

  # Create an option parser
  parser = OptionParser(usage)

  # Add command line options
  parser.add_option('-o', '--output-file',
                    dest='output_file', type="string", default="SPOT.XDS",
                    help='Destination filename for reflections')

  # Parse the arguments
  options, args = parser.parse_args()

  # Print help if no arguments specified, otherwise call function
  if len(args) < 1:
    parser.print_help()
  else:
    # Initialise the script runner
    runner = ScriptRunner(
        reflections_filename=args[0],
        output_filename=options.output_file)

    # Run the script
    runner()
