#!/usr/bin/env python
#
# dials.command_line.index_symops.py
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
from dials.util.options import OptionParser
from dials.util.options import flatten_reflections, flatten_experiments

help_message = '''

This program can be used to (i) generate the symmetry from the experiment
and apply to the input reflections, and (ii) for each symop for that symmetry
attempt to calculate the CC on that operation within the strong spot list.

  dials.index_symops experiment.json indexed.pickle [d_min=3.0]

'''

phil_scope = iotbx.phil.parse("""
d_min = 0
  .type = float
  .help = "Resolution limit to use for analysis"
""", process_includes=True)


def run(args):
  import libtbx.load_env
  from libtbx.utils import Sorry
  usage = "%s [options] experiment.json indexed.pickle" % \
    libtbx.env.dispatcher_name

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
  if len(reflections) == 0 or len(experiments) == 0:
    parser.print_help()
    return
  assert(len(reflections) == 1)
  assert(len(experiments) == 1)
  experiment = experiments[0]
  crystal = experiment.crystal
  reflections = reflections[0]

  from copy import deepcopy

  original_miller_indices = reflections['miller_index']

  space_group = crystal.get_space_group()
  unit_cell = crystal.get_unit_cell()
  from cctbx.crystal import symmetry as crystal_symmetry
  cs = crystal_symmetry(unit_cell, space_group.type().lookup_symbol())

  from cctbx.miller import set as miller_set

  ms = miller_set(cs, original_miller_indices)
  ms = ms.array(reflections['intensity.sum.value'])

  if params.d_min:
    ms = ms.resolution_filter(d_min=params.d_min)

  for smx in space_group.smx():
    reindexed = deepcopy(reflections)
    miller_indices = reflections['miller_index']
    reindexed_miller_indices = sgtbx.change_of_basis_op(smx).apply(
      miller_indices)
    rms = miller_set(cs, reindexed_miller_indices)
    rms = rms.array(reflections['intensity.sum.value'])
    if params.d_min:
      rms = rms.resolution_filter(d_min=params.d_min)
    intensity, intensity_rdx = rms.common_sets(ms)
    cc = intensity.correlation(intensity_rdx).coefficient()

    print '%10s %.3f' % (smx, cc)

if __name__ == '__main__':
  import sys
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run(sys.argv[1:])
