#!/usr/bin/env python
#
# dials.find_spots.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

help_message = '''

This program tries to find strong spots on a sequence of images. The program can
be called with either a "datablock.json" file or a sequence of image files (see
help for dials.import for more information about how images are imported). Spot
finding will be done against each logically grouped set of images given. Strong
pixels will be found on each image and spots will be formed from connected
components. In the case of rotation images, connected component labelling will
be done in 3D.

Once a set of spots have been found, their centroids and intensities will be
calculated. They will then be filtered according to the particular preferences
of the user. The output will be a file (strong.pickle) containing a list of spot
centroids and intensities which can be used in the dials.index program. To view
a list of parameters for spot finding use the --show-config option.

Examples:

  dials.find_spots image1.cbf

  dials.find_spots imager_00*.cbf

  dials.find_spots datablock.json

  dials.find_spots datablock.json output=strong.pickle

'''


class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    from libtbx.phil import parse

    # Set the phil scope
    phil_scope = parse('''

      output = 'strong.pickle'
        .type = str
        .help = "The output filename"

      save_shoeboxes = True
        .type = bool
        .help = "Save the raw pixel values inside the reflection shoeboxes."

      include scope dials.algorithms.peak_finding.spotfinder_factory.phil_scope

    ''', process_includes=True)

    # The script usage
    usage = "usage: %prog [options] [param.phil] "\
            "{datablock.json | image1.file [image2.file ...]}"

    # Initialise the base class
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      epilog=help_message)

  def run(self):
    '''Execute the script.'''
    from dials.util.command_line import Command
    from dials.util.command_line import Importer
    from dials.array_family import flex
    from time import time
    start_time = time()

    # Parse the command line
    params, options, args = self.parser.parse_args(show_diff_phil=True)

    # Try importing the command line arguments
    importer = Importer(args, include=['images', 'datablocks'])

    # Check the unhandled arguments
    if len(importer.unhandled_arguments) > 0:
      print '-' * 80
      print 'The following command line arguments weren\'t handled'
      for arg in importer.unhandled_arguments:
        print '  ' + arg
      exit(1)

    # Ensure we have a data block
    if not importer.datablocks:
      self.config().print_help()
      exit(1)

    if len(importer.datablocks) != 1:
      raise RuntimeError('only 1 datablock can be processed at a time')

    # Loop through all the imagesets and find the strong spots
    reflections = flex.reflection_table.from_observations(
      importer.datablocks[0], params)

    # Delete the shoeboxes
    if not params.save_shoeboxes:
      del reflections['shoebox']

    # Save the reflections to file
    print '\n' + '-' * 80
    Command.start('Saving {0} reflections to {1}'.format(
        len(reflections), params.output))
    reflections.as_pickle(params.output)
    Command.end('Saved {0} reflections to {1}'.format(
        len(reflections), params.output))

    # Print the time
    print time() - start_time


if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
