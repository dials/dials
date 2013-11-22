#!/usr/bin/env python
#
# dials.match_observed_with_predicted.py
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
        type = 'string', default = 'matched.pickle',
        help = 'Set the filename for matched reflections.')

  def main(self, params, options, args):
    '''Execute the script.'''
    from dials.algorithms.integration import ReflectionPredictor
    from dials.algorithms.peak_finding import SpotMatcher
    from dials.util.command_line import Importer, Command
    from dials.model.serialize import dump
    from dials.model.data import ReflectionList

    # Try importing the command line arguments
    importer = Importer(args)
    if len(importer.imagesets) == 0 or len(importer.crystals) == 0:
      self.config().print_help()
      return
    elif len(importer.reflections) == 0:
      self.config().print_help()
      return
    elif len(importer.imagesets) > 1:
      raise RuntimeError("Only one imageset can be processed at a time")
    sweep = importer.imagesets[0]
    crystal = importer.crystals[0]

    # Load the reflections
    Command.start('Loading reflections')
    observed = importer.reflections
    Command.end('Loaded {0} reflections'.format(len(observed)))

    # Create the reflection predictor
    self.predict = ReflectionPredictor()

    # Create the spot matcher
    self.match = SpotMatcher(max_separation=2)

    # Predict the reflections
    predicted = self.predict(sweep, crystal)

    # Match the predictions with the strong spots
    oindex, pindex = self.match(observed, predicted)

    # Copy all of the reflection data for the matched reflections
    Command.start('Creating matches')
    matched = ReflectionList()
    for io, ip in zip(oindex, pindex):
      o = observed[io]
      p = predicted[ip]
      o.miller_index = p.miller_index
      o.rotation_angle = p.rotation_angle
      o.beam_vector = p.beam_vector
      o.image_coord_px = p.image_coord_px
      o.image_coord_mm = p.image_coord_mm
      o.panel_number = p.panel_number
      o.frame_number = p.frame_number
      matched.append(o)
    Command.end('Created {0} matches'.format(len(matched)))

    # Save the reflections to file
    Command.start('Saving {0} reflections to {1}'.format(
        len(matched), options.output_filename))
    dump.reflections(matched, options.output_filename)
    Command.end('Saved {0} reflections to {1}'.format(
        len(matched), options.output_filename))


if __name__ == '__main__':
  script = Script()
  script.run()
