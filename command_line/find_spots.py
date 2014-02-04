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

class Script(ScriptRunner):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''

    # The script usage
    usage = "usage: %prog [options] [param.phil] "\
            "{sweep.json | image1.file [image2.file ...]}"

    # Initialise the base class
    ScriptRunner.__init__(self, usage=usage)

    # Output filename option
    self.config().add_option(
        '-o', '--output-filename',
        dest = 'output_filename',
        type = 'string', default = 'strong.pickle',
        help = 'Set the filename for found strong spots.')

  def main(self, params, options, args):
    '''Execute the script.'''
    from dials.algorithms.peak_finding.spotfinder_factory \
        import SpotFinderFactory
    from dials.algorithms import shoebox
    from dials.model.serialize import load, dump
    from dials.util.command_line import Command
    from dxtbx.imageset import ImageSetFactory
    from dials.util.command_line import Importer
    from dials.array_family import flex

    # Try importing the command line arguments
    importer = Importer(args, include=['images', 'datablocks'])

    # Check the unhandled arguments
    if len(importer.unhandled_arguments) > 0:
      print '-' * 80
      print 'The following command line arguments weren\'t handled'
      for arg in importer.unhandled_arguments:
        print '  ' + arg

    # Ensure we have a data block
    if importer.datablocks is None or len(importer.datablocks) == 0:
      print 'Error: no datablock specified'
      return
    elif len(importer.datablocks) != 1:
      print 'Error: only 1 datablock can be processed at a time'
      return

    # Check the number of imagesets
    imagesets = importer.datablocks[0].extract_imagesets()
    if len(imagesets) == 0:
      print 'Error: no imageset specified'
      return
    elif len(imagesets) != 1:
      print 'Error: only 1 imageset can be processed at a time'
      return
    imageset = imagesets[0]

    # Get the integrator from the input parameters
    print 'Configuring spot finder from input parameters'
    find_spots = SpotFinderFactory.from_parameters(params)

    # Find the strong spots in the sweep
    print 'Finding strong spots'
    reflections = find_spots(imageset)

    # Add the experiment ID
    index = flex.size_t(reflections.nrows(), 0)
    reflections['id'] = index

    # Dump the shoeboxes
    if not params.spotfinder.save_shoeboxes:
      del reflections['shoebox']

    # Save the reflections to file
    Command.start('Saving {0} reflections to {1}'.format(
        len(reflections), options.output_filename))
    dump.reflections(reflections, options.output_filename)
    Command.end('Saved {0} reflections to {1}'.format(
        len(reflections), options.output_filename))


if __name__ == '__main__':
  from time import time
  start_time = time()

  script = Script()
  script.run()

  print "Total time: ", time() - start_time
