#!/usr/bin/env python
#
# dials.model_background.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

from __future__ import division

help_message = '''


'''

# Set the phil scope
from libtbx.phil import parse
phil_scope = parse('''

  output {
    model = 'background.pickle'
      .type = str
      .help = "The output filename"

    log = 'dials.model_background.log'
      .type = str
      .help = "The log filename"

    debug_log = 'dials.model_background.debug.log'
      .type = str
      .help = "The debug log filename"
  }

  verbosity = 1
    .type = int(value_min=0)
    .help = "The verbosity level"

  include scope dials.algorithms.integration.integrator.phil_scope
  include scope dials.algorithms.spot_prediction.reflection_predictor.phil_scope

''', process_includes=True)


class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    import libtbx.load_env

    # The script usage
    usage = "usage: %s [options] [param.phil] "\
            "experiments.json" \
            % libtbx.env.dispatcher_name

    # Initialise the base class
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      epilog=help_message,
      read_experiments=True)

  def run(self):
    '''Execute the script.'''
    from dials.util.command_line import heading
    from dials.array_family import flex
    from dials.util.options import flatten_experiments
    from time import time
    from dials.util import log
    from logging import info, debug
    from libtbx.utils import Sorry
    from dials.algorithms.background.modeller import BackgroundModeller
    start_time = time()

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=False)

    # Configure the logging
    log.config(
      params.verbosity,
      info=params.output.log,
      debug=params.output.debug_log)

    from dials.util.version import dials_version
    info(dials_version())

    # Log the diff phil
    diff_phil = self.parser.diff_phil.as_str()
    if diff_phil is not '':
      info('The following parameters have been modified:\n')
      info(diff_phil)

    # Ensure we have a data block
    experiments = flatten_experiments(params.input.experiments)
    if len(experiments) == 0:
      self.parser.print_help()
      return

    # Only handle a single imageset at once
    imagesets = set(expr.imageset for expr in experiments)
    if len(imagesets) != 1:
      raise Sorry("Can only process a single imageset at a time")

    # Predict the reflections
    info("")
    info("=" * 80)
    info("")
    info(heading("Predicting reflections"))
    info("")
    predicted = flex.reflection_table.from_predictions_multi(
      experiments,
      dmin=params.prediction.d_min,
      dmax=params.prediction.d_max,
      margin=params.prediction.margin,
      force_static=params.prediction.force_static)

    # Create the modeller
    modeller = BackgroundModeller(experiments, predicted, params)
    model = modeller.compute()

    # Save the background model
    info("Saving background model to %s" % params.output.model)
    with open(params.output.model, "w") as outfile:
      import cPickle as pickle
      pickle.dump(model, outfile, protocol=pickle.HIGHEST_PROTOCOL)

    # Print the time
    info("Time Taken: %f" % (time() - start_time))


if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
