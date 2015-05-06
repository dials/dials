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

help_message = '''

This program can be used to re-index an experiments.json and/or indexed.pickle
file from one setting to another. The change of basis operator can be
provided in h,k,l, or a,b,c or x,y,z conventions. By default the change of
basis operator will also be applied to the space group in the experiments.json
file, however, optionally, a space group (including setting) to be applied
AFTER applying the change of basis operator can be provided.

Examples::

  dials.reindex experiments.json change_of_basis_op=b+c,a+c,a+b

  dials.reindex indexed.pickle change_of_basis_op=-b,a+b+2*c,-a

  dials.reindex experiments.json index.pickle change_of_basis_op=l,h,k

'''

phil_scope = iotbx.phil.parse("""
change_of_basis_op = a,b,c
  .type = str
space_group = None
  .type = space_group
  .help = "The space group to be applied AFTER applying the change of basis "
           "operator."

output {
  experiments = reindexed_experiments.json
    .type = str
    .help = "The filename for reindexed experimental models"

  reflections = reindexed_reflections.pickle
    .type = str
    .help = "The filename for reindexed reflections"
}
""", process_includes=True)


def run(args):
  import libtbx.load_env
  from libtbx.utils import Sorry
  usage = "%s [options] experiments.json indexed.pickle" %libtbx.env.dispatcher_name

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_reflections=True,
    read_experiments=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)

  reflections = flatten_reflections(params.input.reflections)
  experiments = flatten_experiments(params.input.experiments)
  if len(experiments) == 0 and len(reflections) == 0:
    parser.print_help()
    return
  elif len(experiments) > 1:
    raise Sorry("Only one Experiment can be processed at a time")
  if params.change_of_basis_op is None:
    raise Sorry("Please provide a change_of_basis_op.")

  change_of_basis_op = sgtbx.change_of_basis_op(params.change_of_basis_op)

  if len(experiments):
    experiment = experiments[0]
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

    dump.experiment_list(experiments, params.output.experiments)

  if len(reflections):
    assert(len(reflections) == 1)
    reflections = reflections[0]

    miller_indices = reflections['miller_index']
    miller_indices_reindexed = change_of_basis_op.apply(miller_indices)
    reflections['miller_index'] = miller_indices_reindexed

    easy_pickle.dump(params.output.reflections, reflections)


if __name__ == '__main__':
  import sys
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run(sys.argv[1:])
