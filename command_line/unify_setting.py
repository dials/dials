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

from __future__ import absolute_import, division, print_function

import iotbx.phil
from dials.util.options import OptionParser
from libtbx.phil import command_line

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
  from scitbx.math import r3_rotation_axis_and_angle_from_matrix

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

  for j, experiment in enumerate(experiments):
    goniometer = experiment.data.goniometers()[0]
    F = matrix.sqr(goniometer.get_fixed_rotation())
    crystal = experiment.data.crystals()[0]
    U = matrix.sqr(crystal.get_U())
    B = matrix.sqr(crystal.get_B())
    UB = F * U * B
    UBt = UB.transpose().elems
    a, b, c = matrix.col(UBt[0:3]), matrix.col(UBt[3:6]), matrix.col(UBt[6:9])
    axis = matrix.col(goniometer.get_rotation_axis())
    from math import pi
    r2d = 180 / pi
    abc = [a, b, c]
    abc_names = 'abc'
    distances = [(r2d * (min(axis.angle(_a), pi - axis.angle(_a))), k)
                 for k, _a in enumerate(abc)]
    close = sorted(distances)[0]
    if reference_U is None:
      reference_U = U
      reference_space_group = lattice_symmetry_group(crystal.get_unit_cell(),
                                                     max_delta=0.0)
      print('%s possible lattice ops' % len(reference_space_group.all_ops()))

    print('Experiment %d' % j)
    print('Closest (original) axis: %s* %.2f' % \
      (abc_names[close[1]], close[0]))

    results = []
    for op in reference_space_group.all_ops():
      R = B * matrix.sqr(op.r().as_double()).transpose() * B.inverse()
      relative = (U * R).inverse() * reference_U
      rot = r3_rotation_axis_and_angle_from_matrix(relative)
      results.append((abs(rot.angle()), op.r().as_hkl(), rot))
    results.sort()
    print('Best reindex op for experiment %d: %12s (%.3f)' % \
      (j, results[0][1], 180.0 * results[0][2].angle() / pi))

    if results[0][0] > (5 * pi / 180.0):
      print('Rotation: axis: %.4f %.4f %.4f' % results[0][2].axis)
      print('          angle: %.4f degrees' % \
        (180.0 * results[0][2].angle() / pi))

if __name__ == '__main__':
  import sys
  from libtbx.utils import show_times_at_exit
  run(sys.argv[1:])
