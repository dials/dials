#!/usr/bin/env python
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# dials.command_line.reindex.py
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Richard Gildea
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

from libtbx.phil import command_line
from libtbx import easy_pickle
import iotbx.phil
from cctbx import sgtbx
from dxtbx.model.crystal import crystal_model
from dxtbx.serialize import dump
# from dials.util.command_line import Importer
from dials.util.options import OptionParser
from dials.util.options import flatten_reflections, flatten_experiments

master_phil_scope = iotbx.phil.parse("""
change_of_basis_op = None
  .type = str
  .optional = False
space_group = None
  .type = space_group
  .help = "The space group to be applied AFTER applying the change of basis "
           "operator."
""", process_includes=True)

master_params = master_phil_scope.fetch().extract()


def run(args):

  parser = OptionParser(
    phil=master_phil_scope,
    read_reflections=True,
    read_experiments=True,
    check_format=False)

  params, options = parser.parse_args(show_diff_phil=True)

  reflections = flatten_reflections(params.input.reflections)
  experiments = flatten_experiments(params.input.experiments)
  if len(experiments) == 0:
    print "No ExperimentList could be constructed"
    return
  elif len(experiments) > 1:
    raise RuntimeError("Only one Experiment can be processed at a time")

  # importer = Importer(args, check_format=False)
  # if len(importer.experiments) == 0:
  #   print "No ExperimentList could be constructed"
  #   return
  # elif len(importer.experiments) > 1:
  #   raise RuntimeError("Only one Experiment can be processed at a time")
  # experiments = importer.experiments
  experiment = experiments[0]
  assert(len(reflections) == 1)
  reflections = reflections[0]
  # assert len(importer.reflections) == 1
  # reflections = importer.reflections[0]
  # args = importer.unhandled_arguments
  # cmd_line = command_line.argument_interpreter(master_params=master_phil_scope)
  # working_phil = cmd_line.process_and_fetch(args=args)
  # working_phil.show()
  # params = working_phil.extract()
  parser.phil.show()
  assert params.change_of_basis_op is not None

  change_of_basis_op = sgtbx.change_of_basis_op(params.change_of_basis_op)
  cryst_orig = experiment.crystal
  cryst_reindexed = cryst_orig.change_basis(change_of_basis_op)
  if params.space_group is not None:
    a, b, c = cryst_reindexed.get_real_space_vectors()
    cryst_reindexed = crystal_model(
      a, b, c, space_group=params.space_group.group())
  experiment.crystal = cryst_reindexed

  print "Old crystal:"
  print cryst_orig
  print
  print "New crystal:"
  print cryst_reindexed
  print

  miller_indices = reflections['miller_index']
  miller_indices_reindexed = change_of_basis_op.apply(miller_indices)
  reflections['miller_index'] = miller_indices_reindexed

  easy_pickle.dump('reflections_reindexed.pickle', reflections)
  dump.experiment_list(experiments, 'experiments_reindexed.json')


if __name__ == '__main__':
  import sys
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run(sys.argv[1:])
