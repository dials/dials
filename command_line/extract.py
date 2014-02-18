#!/usr/bin/env python
#
# dials.extract.py
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

    # The block length
    self.config().add_option(
        '-n', '--num-blocks',
        dest = 'num_blocks',
        type = 'int', default = 1,
        help = 'Set the number of blocks')

    # Output filename option
    self.config().add_option(
        '-o', '--output-filename',
        dest = 'output_filename',
        type = 'string', default = 'extracted.tar',
        help = 'Set the filename for the extracted spots.')

    # Output filename option
    self.config().add_option(
        '--force-static',
        dest = 'force_static',
        action = "store_true", default = False,
        help = 'For a scan varying model force static prediction.')

  def main(self, params, options, args):
    '''Execute the script.'''
    from dials.model.serialize import load, dump
    from dials.util.command_line import Command
    from dials.util.command_line import Importer
    from dials.algorithms.shoebox import ProfileBlockExtractor
    from dials.algorithms.shoebox import BBoxCalculator
    from dials.array_family import flex
    from dials.model.data import ReflectionList

    # Check the unhandled arguments
    importer = Importer(args, include=['experiments'])
    if len(importer.unhandled_arguments) > 0:
      print '-' * 80
      print 'The following command line arguments weren\'t handled'
      for arg in importer.unhandled_arguments:
        print '  ' + arg

    # Check the number of experiments
    if importer.experiments is None or len(importer.experiments) == 0:
      print 'Error: no experiment list specified'
      return
    elif len(importer.experiments) > 1:
      print 'Error: only 1 experiment currently supported'
      return

    # Populate the reflection table with predictions
    predicted = flex.reflection_table.from_predictions(
      importer.experiments,
      force_static=options.force_static)

    # Get the bbox nsigma
    n_sigma = params.shoebox.n_sigma

    # Calculate the bounding boxes
    predicted.compute_bbox(importer.experiments[0], n_sigma)

    # Create the profile block extractor
    extract = ProfileBlockExtractor(
      importer.experiments[0].imageset, predicted,
      options.num_blocks, options.output_filename)


if __name__ == '__main__':
  script = Script()
  script.run()
