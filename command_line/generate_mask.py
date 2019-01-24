#!/usr/bin/env python
#
# generate_mask.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division, print_function

from iotbx.phil import parse

help_message = '''

This program is used to generate mask to specify which pixels should be
considered "invalid" during spot finding and integration. It provides a few
options to create simple masks using the detector trusted range, or from
simple shapes or by setting different resolution ranges.

Examples::

  dials.generate_mask experiments.json border=5

  dials.generate_mask experiments.json \\
    untrusted.rectangle=50,100,50,100 \\
    untrusted.circle=200,200,100

  dials.generate_mask experiments.json resolution.d_max=2.00

'''

phil_scope = parse("""
  output {
    mask = mask.pickle
      .type = path
      .help = "Name of output mask file"
    experiments = None
      .type = path
      .help = "Save the modified experiments. (usually only modified with the"
              "generated pixel mask)"
  }

  include scope dials.util.masking.phil_scope
""", process_includes=True)


class Script(object):
  ''' A class to encapsulate the script. '''

  def __init__(self):
    ''' Initialise the script. '''
    from dials.util.options import OptionParser
    import libtbx.load_env

    # Create the parser
    usage = "usage: %s [options] experiments.json" % libtbx.env.dispatcher_name
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      epilog=help_message,
      read_experiments=True)

  def run(self):
    ''' Run the script. '''
    from dials.util.masking import MaskGenerator
    from dials.util.options import flatten_experiments
    from libtbx.utils import Sorry
    import six.moves.cPickle as pickle
    from dials.util import log
    from dxtbx.format.image import ImageBool

    # Parse the command line arguments
    params, options = self.parser.parse_args(show_diff_phil=True)
    experiments = flatten_experiments(params.input.experiments)

    # Configure logging
    log.config()

    # Check number of args
    if len(experiments) == 0:
      self.parser.print_help()
      return

    if len(experiments) != 1:
      raise Sorry('exactly 1 experiments must be specified')
    imagesets = experiments.imagesets()
    if len(imagesets) > 1:
      # Check beams (for resolution) and detectors are equivalent in each case
      # otherwise the mask may not be appropriate across all imagesets
      detectors = experiments.detectors()
      beams = experiments.beams()
      for d in detectors[1:]:
        if not d.is_similar_to(detectors[0]):
          raise Sorry('multiple imagesets are present, but their detector'
                      ' models differ')
      for b in beams[1:]:
        if not b.is_similar_to(beams[0]):
          raise Sorry('multiple imagesets are present, but their beam'
                      ' models differ')

    imageset = imagesets[0]

    # Generate the mask
    generator = MaskGenerator(params)
    mask = generator.generate(imageset)

    # Save the mask to file
    print("Writing mask to %s" % params.output.mask)
    with open(params.output.mask, "wb") as fh:
      pickle.dump(mask, fh)

    # Save the experiment list
    if params.output.experiments is not None:
      for imageset in imagesets:
        imageset.external_lookup.mask.data = ImageBool(mask)
        imageset.external_lookup.mask.filename = params.output.mask
      from dxtbx.model.experiment_list import ExperimentListDumper
      print('Saving experiments to {0}'.format(
        params.output.experiments))
      dump = ExperimentListDumper(experiments)
      dump.as_file(params.output.experiments)

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
