#!/usr/bin/env python
# LIBTBX_SET_DISPATCHER_NAME dev.dials.unify_setting
#
# dials.command_line.unify_setting.py
#
#  Copyright (C) 2015 Diamond Light Source
#
#  Author: Graeme Winter
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

from libtbx.phil import command_line
import iotbx.phil
from cctbx import sgtbx
from dials.util.options import OptionParser

help_message = '''

This program is a jiffy to work out how to reindex orientation matrices to
match one another, when multiple fixed_rotations are in play.

'''

phil_scope = iotbx.phil.parse("""
""", process_includes=True)

def run(args):
  import libtbx.load_env
  from scitbx import matrix
  from cctbx.sgtbx import lattice_symmetry_group
  usage = "%s [options] experiment_0.json ..." % \
    libtbx.env.dispatcher_name

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_experiments=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)
  experiments = params.input.experiments

  # check input
  space_group = None
  for experiment in experiments:
    assert(len(experiment.data.goniometers()) == 1)
    assert(len(experiment.data.crystals()) == 1)
    crystal = experiment.data.crystals()[0]
    if space_group is None:
      space_group = crystal.get_space_group()
    else:
      assert(crystal.get_space_group() == space_group)

  reference_U = None
  reference_space_group = None
  i3 = matrix.sqr((1, 0, 0, 0, 1, 0, 0, 0, 1))

  for j, experiment in enumerate(experiments):
    goniometer = experiment.data.goniometers()[0]
    crystal = experiment.data.crystals()[0]
    fixed = matrix.sqr(goniometer.get_fixed_rotation())
    U = matrix.sqr(crystal.get_U())
    B = matrix.sqr(crystal.get_B())
    if reference_U is None:
      reference_U = U
      reference_space_group = lattice_symmetry_group(crystal.get_unit_cell(),
                                                     max_delta=0.0)
      print '%s possible lattice symops' % len(reference_space_group.all_ops())

    results = []
    for op in reference_space_group.all_ops():
      R = B * matrix.sqr(op.r().as_double()).transpose() * B.inverse()
      nearly_i3 = (U * R).inverse() * reference_U
      score = sum([abs(_n - _i) for (_n, _i) in zip(nearly_i3, i3)])
      results.append((score, op.r().as_hkl()))
    results.sort()
    print 'Best reindex op for experiment %d: %12s (%.3f)' % \
      (j, results[0][1], results[0][0])






if __name__ == '__main__':
  import sys
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run(sys.argv[1:])
