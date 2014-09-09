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

class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    from libtbx.phil import parse

    # The script usage
    usage = "usage: %prog [options] [param.phil] "\
            "{sweep.json | image1.file [image2.file ...]}"

    phil_scope = parse('''
      output = predicted.pickle
        .type = str
        .help = "The filename for the predicted reflections"

      force_static = False
        .type = bool
        .help = "For a scan varying model, force static prediction"

      buffer_size = 0
        .type = int
        .help = "Calculate predictions within a buffer zone of n images either"
                "size of the scan"

      dmin = None
        .type = float
        .help = "Minimum d-spacing of predicted reflections"
    ''')

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=self.phil_scope())

  def run(self):
    '''Execute the script.'''
    from dials.model.serialize import load, dump
    from dials.util.command_line import Command
    from dials.util.command_line import Importer
    from dials.array_family import flex
    from dials.framework.registry import Registry
    from dials.algorithms.profile_model.profile_model import ProfileModel

    # Parse the command line
    params, options, args = self.parser.parse_args()

    # Check the unhandled arguments
    importer = Importer(args, include=['experiments'], check_format=False)
    if len(importer.unhandled_arguments) > 0:
      print '-' * 80
      print 'The following command line arguments weren\'t handled'
      for arg in importer.unhandled_arguments:
        print '  ' + arg
      exit(1)

    # Check the number of experiments
    if importer.experiments is None or len(importer.experiments) == 0:
      self.config().print_help()
      return

    predicted_all = flex.reflection_table()

    for i_expt, expt in enumerate(importer.experiments):
      if params.buffer_size > 0:
        # Hack to make the predicter predict reflections outside of the range
        # of the scan
        scan = expt.scan
        image_range = scan.get_image_range()
        oscillation = scan.get_oscillation()
        scan.set_image_range((image_range[0]-params.buffer_size,
                              image_range[1]+params.buffer_size))
        scan.set_oscillation((oscillation[0]-params.buffer_size*oscillation[1],
                              oscillation[1]))

      # Populate the reflection table with predictions
      predicted = flex.reflection_table.from_predictions(
        expt,
        force_static=params.force_static,
        dmin=params.dmin)
      predicted['id'] = flex.int(len(predicted), i_expt)

      # Compute the bounding box
      registry = Registry()
      n_sigma = registry.params().integration.shoebox.n_sigma
      sigma_b = registry.params().integration.shoebox.sigma_b
      sigma_m = registry.params().integration.shoebox.sigma_m
      if sigma_b is not None and sigma_m is not None:
        import math
        d2r = math.pi / 180.0
        profile_model = ProfileModel(n_sigma, sigma_b*d2r, sigma_m*d2r)
        predicted.compute_bbox(expt, profile_model)
      predicted_all.extend(predicted)

    # Save the reflections to file
    Command.start('Saving {0} reflections to {1}'.format(
        len(predicted_all), params.output))
    predicted_all.as_pickle(params.output)
    Command.end('Saved {0} reflections to {1}'.format(
        len(predicted_all), params.output))


if __name__ == '__main__':
  script = Script()
  script.run()
