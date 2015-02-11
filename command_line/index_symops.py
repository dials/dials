#!/usr/bin/env python
#
# dials.command_line.index_symops.py
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
from dials.util.options import flatten_reflections, flatten_experiments

help_message = '''

This program can be used to (i) generate the symmetry from the experiment
and apply to the input reflections, and (ii) for each symop for that symmetry
attempt to calculate the CC on that operation within the strong spot list.

  dials.index_symops experiment.json indexed.pickle [d_min=3.0] [d_max=10.0]

'''

phil_scope = iotbx.phil.parse("""
d_min = 0
  .type = float
  .help = "High resolution limit to use for analysis"
d_max = 0
  .type = float
  .help = "Low resolution limit to use for analysis"
symop_threshold = 0
  .type = float
  .help = "Thereshold above which we consider a symmetry operator true."
""", process_includes=True)


def run(args):
  import libtbx.load_env
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
  from dials.array_family import flex

  original_miller_indices = reflections['miller_index']

  space_group = crystal.get_space_group()
  unit_cell = crystal.get_unit_cell()
  from cctbx.crystal import symmetry as crystal_symmetry
  cs = crystal_symmetry(unit_cell, space_group.type().lookup_symbol())

  from cctbx.miller import set as miller_set

  ms = miller_set(cs, original_miller_indices)
  ms = ms.array(reflections['intensity.sum.value'] /
                flex.sqrt(reflections['intensity.sum.variance']))

  if params.d_min or params.d_max:
    ms = ms.resolution_filter(d_min=params.d_min, d_max=params.d_max)

  print '%d reflections analysed' % ms.size()
  print '%10s %6s %5s' % ('Symop', 'Nref', 'CC')

  true_symops = []

  for smx in space_group.smx():
    reindexed = deepcopy(reflections)
    miller_indices = reflections['miller_index']
    reindexed_miller_indices = sgtbx.change_of_basis_op(smx).apply(
      miller_indices)
    rms = miller_set(cs, reindexed_miller_indices)
    rms = rms.array(reflections['intensity.sum.value'] /
                    flex.sqrt(reflections['intensity.sum.variance']))
    if params.d_min or params.d_max:
      rms = rms.resolution_filter(d_min=params.d_min, d_max=params.d_max)
    intensity, intensity_rdx = rms.common_sets(ms)
    cc = intensity.correlation(intensity_rdx).coefficient()

    accept = ''

    if params.symop_threshold:
      if cc > params.symop_threshold:
        true_symops.append(smx)
        accept = '***'

    print '%10s %6d %.3f %s' % (smx, intensity.size(), cc, accept)

  if params.symop_threshold:
    from cctbx.sgtbx import space_group
    sg = space_group()
    for symop in true_symops:
      sg = sg.expand_smx(symop)
    sg_symbols = sg.match_tabulated_settings()
    print 'Derived point group from symmetry operations: %s' % \
      sg_symbols.hermann_mauguin()

if __name__ == '__main__':
  import sys
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run(sys.argv[1:])
