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
from dials.model.serialize import dump
from dials.util.command_line import Importer

master_phil_scope = iotbx.phil.parse("""
change_of_basis_op = None
  .type = str
  .optional = False
""", process_includes=True)

master_params = master_phil_scope.fetch().extract()


def run(args):
  importer = Importer(args, check_format=False)
  if len(importer.experiments) == 0:
    print "No ExperimentList could be constructed"
    return
  elif len(importer.experiments) > 1:
    raise RuntimeError("Only one Experiment can be processed at a time")
  experiments = importer.experiments
  experiment = experiments[0]
  assert len(importer.reflections) == 1
  reflections = importer.reflections[0]
  args = importer.unhandled_arguments
  cmd_line = command_line.argument_interpreter(master_params=master_phil_scope)
  working_phil = cmd_line.process_and_fetch(args=args)
  working_phil.show()
  params = working_phil.extract()
  assert params.change_of_basis_op is not None

  change_of_basis_op = sgtbx.change_of_basis_op(params.change_of_basis_op)
  crystal_model = experiment.crystal
  crystal_model_reindexed = crystal_model.change_basis(change_of_basis_op)
  experiment.crystal = crystal_model_reindexed

  print "Old crystal:"
  print crystal_model
  print
  print "New crystal:"
  print crystal_model_reindexed
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
