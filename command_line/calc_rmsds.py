#!/usr/bin/env python
#
# calc_rmsds.py
#
#  Copyright (C) 2016 Diamond Light Source
#
#  Author: Richard Gildea
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

# DIALS_ENABLE_COMMAND_LINE_COMPLETION
from __future__ import absolute_import, division, print_function

import logging

import libtbx.load_env

logger = logging.getLogger('dials.command_line.calc_rmsds')

help_message = '''

'''

# The phil scope
from libtbx.phil import parse
phil_scope = parse('''
''', process_includes=True)

from dials.command_line.refine import phil_scope


def run(args):
  from dials.util import log
  usage = "%s experiments.json indexed.pickle [options]" %libtbx.env.dispatcher_name


  from dials.util.options import OptionParser
  from dials.util.options import flatten_reflections
  from dials.util.options import flatten_experiments

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_experiments=True,
    read_reflections=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=False)

  # Configure the logging
  #log.config(info=params.output.log, debug=params.output.debug_log)

  from dials.util.version import dials_version
  logger.info(dials_version())

  # Log the diff phil
  diff_phil = parser.diff_phil.as_str()
  if diff_phil is not '':
    logger.info('The following parameters have been modified:\n')
    logger.info(diff_phil)

  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)
  assert(len(reflections) == 1)
  reflections = reflections[0]

  if len(experiments) == 0:
    parser.print_help()
    return

  #from dials.command_line import refine
  #params = refine.phil_scope.extract()
  indexed_reflections = reflections.select(reflections['id'] > -1)
  from dials.algorithms.refinement import RefinerFactory
  refiner = RefinerFactory.from_parameters_data_experiments(
    params, indexed_reflections, experiments)
  #refiner.run()
  rmsds = refiner.rmsds()
  import math
  xy_rmsds = math.sqrt(rmsds[0]**2 + rmsds[1]**2)

  print(rmsds)



  return

if __name__ == '__main__':
  import sys
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run(sys.argv[1:])
