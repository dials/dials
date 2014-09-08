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
from dials.util.script import ScriptRunner

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

  dials.find_spots datablock.json -o strong.pickle

'''


class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser

    # The script usage
    usage = "usage: %prog [options] [param.phil] "\
            "{datablock.json | image1.file [image2.file ...]}"

    # Initialise the base class
    self.parser = OptionParser(
      usage=usage,
      phil=self.phil_scope(),
      epilog=help_message)

    # Add an option to show configuration parameters
    self.parser.add_option(
      '-c',
      action='count',
      default=0,
      dest='show_config',
      help='Show the configuration parameters.')

  def phil_scope(self):
    ''' Get the phil scope. '''
    from libtbx.phil import parse
    new_phil_scope = parse('''
      spotfinder {

        output = 'strong.pickle'
          .type = str
          .help = "The output filename"

        save_shoeboxes = True
          .type = bool
          .help = "Save the raw pixel values inside the reflection shoeboxes."
      }

      include scope dials.algorithms.peak_finding.spotfinder_factory.phil_scope

    ''', process_includes=True)
    return new_phil_scope

  def run(self):
    '''Execute the script.'''
    from dials.util.command_line import Command
    from dials.util.command_line import Importer
    from dials.array_family import flex

    # Parse the command line
    params, options, args = self.parser.parse_args()

    # Show config
    if options.show_config > 0:
      self.parser.print_phil(attributes_level=options.show_config-1)
      return

    # Check the number of command line arguments
    if len(args) != 1:
      self.parser.print_help()
      return

    # Try importing the command line arguments
    importer = Importer(args, include=['images', 'datablocks'])

    # Print the diff phil
    diff_phil_str = self.parser.diff_phil().as_str()
    print 'Finding spots with the following user specified parameters:\n'
    if (diff_phil_str is not ''):
      print diff_phil_str
      print ''
    else:
      print 'All parameters set to defaults\n'

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
    if not params.spotfinder.save_shoeboxes:
      del reflections['shoebox']

    # Save the reflections to file
    print '\n' + '-' * 80
    Command.start('Saving {0} reflections to {1}'.format(
        len(reflections), params.spotfinder.output))
    reflections.as_pickle(params.spotfinder.output)
    Command.end('Saved {0} reflections to {1}'.format(
        len(reflections), params.spotfinder.output))


if __name__ == '__main__':
  from time import time
  start_time = time()

  script = Script()
  script.run()

  print "Total time: ", time() - start_time
