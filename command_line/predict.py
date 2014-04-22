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

    # Populate the reflection table with predictions
    predicted = flex.reflection_table.from_predictions(
      importer.experiments[0],
      force_static=options.force_static)
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
