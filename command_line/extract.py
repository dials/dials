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
      output = shoebox.dat
        .type = str
        .help = "The filename for the shoebox output"

      force_static = False
        .type = bool
        .help = "For a scan varying model, force static prediction"

      dmin = None
        .type = float
        .help = "Minimum d-spacing of predicted reflections"

        include scope dials.algorithms.profile_model.profile_model.phil_scope
    ''', process_includes=True)

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope)

  def run(self):
    '''Execute the script.'''
    from dials.model.serialize import load, dump
    from dials.util.command_line import Command
    from dials.util.command_line import Importer
    from dials.array_family import flex
    from dials.model.serialize import extract_shoeboxes_to_file
    from dials.algorithms.profile_model.profile_model import ProfileModelList
    from math import pi

    # Parse the command line
    params, options, args = self.parser.parse_args(show_diff_phil=True)

    # Check the unhandled arguments
    importer = Importer(args, include=['experiments'])
    if len(importer.unhandled_arguments) > 0:
      print '-' * 80
      print 'The following command line arguments weren\'t handled'
      for arg in importer.unhandled_arguments:
        print '  ' + arg

    # Check the number of experiments
    if importer.experiments is None or len(importer.experiments) == 0:
      self.config().print_help()
      return
    elif len(importer.experiments) > 1:
      print 'Error: only 1 experiment currently supported'
      return

    # Populate the reflection table with predictions
    predicted = flex.reflection_table.from_predictions(
      importer.experiments[0],
      force_static=params.force_static,
      dmin=params.dmin)
    predicted['id'] = flex.size_t(len(predicted), 0)

    # Get the bbox nsigma
    profile_model = ProfileModelList.load(params)

    # Calculate the bounding boxes
    predicted.compute_bbox(importer.experiments, profile_model)

    # TODO Need to save out reflections
    z = predicted['xyzcal.px'].parts()[2]
    index = sorted(range(len(z)), key=lambda x: z[x])
    predicted.reorder(flex.size_t(index))

    # Extract the shoeboxes to file
    extract_shoeboxes_to_file(
      params.output,
      importer.experiments[0].imageset,
      predicted)


if __name__ == '__main__':
  script = Script()
  script.run()
