#!/usr/bin/env python
#
# dials.predict.py
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
        type = 'string', default = 'predicted.pickle',
        help = 'Set the filename for the predicted spots.')

  def main(self, params, options, args):
    '''Execute the script.'''
    from dials.model.serialize import load, dump
    from dials.util.command_line import Command
    from dials.util.command_line import Importer
    from dials.algorithms.integration import ReflectionPredictor

    # Try importing the command line arguments
    importer = Importer(args)
    if len(importer.imagesets) == 0 or len(importer.imagesets) == 0:
      self.config().print_help()
      return
    elif len(importer.imagesets) > 1:
      raise RuntimeError("Only one imageset can be processed at a time")

    sweep = importer.imagesets[0]
    crystal = importer.crystals[0]

    # Predict the reflections
    predict = ReflectionPredictor()
    predicted = predict(sweep, crystal)

    # Save the reflections to file
    Command.start('Saving {0} reflections to {1}'.format(
        len(predicted), options.output_filename))
    dump.reflections(predicted, options.output_filename)
    Command.end('Saved {0} reflections to {1}'.format(
        len(predicted), options.output_filename))


if __name__ == '__main__':
  script = Script()
  script.run()
