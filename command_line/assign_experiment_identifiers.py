#!/usr/bin/env python
# LIBTBX_SET_DISPATCHER_NAME dev.dials.assign_experiment_identifiers
# coding: utf-8
from __future__ import absolute_import, division, print_function

import sys
from libtbx import phil
from dials.util import halraiser
from dials.util.options import OptionParser, flatten_reflections,\
  flatten_experiments
from dials.array_family import flex
from dials.algorithms.scaling.scaling_utilities import assign_unique_identifiers,\
  parse_multiple_datasets, save_experiments, save_reflections

help_message = """Command line script which assigns experiment identifiers
to reflections and experiments and saves them back to disk.
"""

phil_scope = phil.parse('''
  output {
    reflections = assigned_reflections.pickle
      .type = str
    experiments = assigned_experiments.json
      .type = str
  }
''')

if __name__ == "__main__":
  try:
    #Parse the command line and flatten reflections, experiments
    usage = '''Usage: dev.dials.assign_experiment_identifiers reflections.pickle
experiments.json'''
    optionparser = OptionParser(usage=usage, read_experiments=True,
      read_reflections=True, read_datablocks=False, phil=phil_scope,
      check_format=False, epilog=help_message)
    params, _ = optionparser.parse_args(show_diff_phil=False)
    if not params.input.experiments or not params.input.reflections:
      optionparser.print_help()
      sys.exit()
    reflections = flatten_reflections(params.input.reflections)
    experiments = flatten_experiments(params.input.experiments)
    reflections = parse_multiple_datasets(reflections)
    experiments, reflections = assign_unique_identifiers(experiments, reflections)
    save_experiments(experiments, params.output.experiments)
    joint_table = flex.reflection_table()
    for reflection_table in reflections:
      joint_table.extend(reflection_table)
    save_reflections(joint_table, params.output.reflections)

  except Exception as e:
    halraiser(e)
