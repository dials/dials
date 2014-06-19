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

    # Output filename option
    self.config().add_option(
        '--force-static',
        dest = 'force_static',
        action = "store_true", default = False,
        help = 'For a scan varying model force static prediction.')

    # buffer size option
    self.config().add_option(
        '-b', '--buffer-size',
        dest = 'buffer_size',
        type = 'int', default = 0,
        help = 'Calculate predictions within a buffer zone of n images either '
               'side of the scan.')

    self.config().add_option(
        '--dmin',
        dest = 'dmin',
        type = 'int', default = None,
        help = 'Minimum d-spacing of predicted reflections.')


  def main(self, params, options, args):
    '''Execute the script.'''
    from dials.model.serialize import load, dump
    from dials.util.command_line import Command
    from dials.util.command_line import Importer
    from dials.array_family import flex
    from dials.framework.registry import Registry

    # Check the unhandled arguments
    importer = Importer(args, include=['experiments'], check_format=False)
    if len(importer.unhandled_arguments) > 0:
      print '-' * 80
      print 'The following command line arguments weren\'t handled'
      for arg in importer.unhandled_arguments:
        print '  ' + arg

    # Check the number of experiments
    if importer.experiments is None or len(importer.experiments) == 0:
      print 'Error: no experiment list specified'
      return
    assert(len(importer.experiments) == 1)

    if options.buffer_size > 0:
      # Hack to make the predicter predict reflections outside of the range
      # of the scan
      scan = importer.experiments[0].scan
      image_range = scan.get_image_range()
      oscillation = scan.get_oscillation()
      scan.set_image_range((image_range[0]-options.buffer_size,
                            image_range[1]+options.buffer_size))
      scan.set_oscillation((oscillation[0]-options.buffer_size*oscillation[1],
                            oscillation[1]))

    # Populate the reflection table with predictions
    predicted = flex.reflection_table.from_predictions(
      importer.experiments[0],
      force_static=options.force_static,
      dmin=options.dmin)
    predicted['id'] = flex.size_t(len(predicted), 0)

    # Compute the bounding box
    registry = Registry()
    n_sigma = registry.params().integration.shoebox.n_sigma
    sigma_b = registry.params().integration.shoebox.sigma_b
    sigma_m = registry.params().integration.shoebox.sigma_m
    if sigma_b is not None and sigma_m is not None:
      import math
      d2r = math.pi / 180.0
      predicted.compute_bbox(
        importer.experiments[-1],
        nsigma=n_sigma,
        sigma_d=sigma_b * d2r,
        sigma_m=sigma_m * d2r)

    # Save the reflections to file
    Command.start('Saving {0} reflections to {1}'.format(
        len(predicted), options.output_filename))
    predicted.as_pickle(options.output_filename)
    Command.end('Saved {0} reflections to {1}'.format(
        len(predicted), options.output_filename))


if __name__ == '__main__':
  script = Script()
  script.run()
