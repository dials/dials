#!/usr/bin/env python
#
# dials.command_line.check_indexing_symmetry.py
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

  dials.check_indexing_symmetry experiment.json indexed.pickle \
    [d_min=3.0] [d_max=10.0]

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
  .help = "Threshold above which we consider a symmetry operator true."
grid_search_scope = 0
  .type = int
  .help = "Search scope for testing misindexing on h, k, l."
asu = False
  .type = bool
  .help = "Perform search comparing within ASU (assumes input symm)"
normalise = False
  .type = bool
  .help = "Normalise intensities before calculating correlation coefficients."
normalise_bins = 0
  .type = int
  .help = "Number of resolution bins for normalisation"
""", process_includes=True)


def get_symop_correlation_coefficients(miller_array):
  from copy import deepcopy
  from scitbx.array_family import flex
  from cctbx import miller
  corr_coeffs = flex.double()
  n_refs = flex.int()
  space_group = miller_array.space_group()
  for smx in space_group.smx():
    miller_indices = miller_array.indices()
    reindexed_miller_indices = sgtbx.change_of_basis_op(smx).apply(
      miller_indices)
    rms = miller.set(miller_array, reindexed_miller_indices)
    rms = rms.array(miller_array.data())
    intensity, intensity_rdx = rms.common_sets(miller_array)
    cc = intensity.correlation(intensity_rdx).coefficient()
    corr_coeffs.append(cc)
    n_refs.append(intensity.size())
  return corr_coeffs, n_refs

def normalise_intensities(miller_array, n_bins=10):
  miller_array.setup_binner(n_bins=n_bins)
  nomalisations = miller_array.amplitude_quasi_normalisations()
  miller_array = miller_array.customized_copy(
    data=miller_array.data()/nomalisations.data())
  return miller_array

def test_crystal_pointgroup_symmetry(reflections, experiment, params):
  crystal = experiment.crystal

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
    d_spacings = ms.d_spacings().data()
    sel = (d_spacings >= params.d_min) & (d_spacings <= params.d_max)
    ms = ms.select(sel)
    reflections = reflections.select(sel)

  if params.normalise:
    if params.normalise_bins:
      ms = normalise_intensities(ms, n_bins=params.normalise_bins)
    else:
      ms = normalise_intensities(ms)
  print 'Check symmetry operations on %d reflections:' % ms.size()
  print ''
  print '%10s %6s %5s' % ('Symop', 'Nref', 'CC')

  true_symops = []

  ccs, n_refs = get_symop_correlation_coefficients(ms)

  for smx, cc, n_ref in zip(space_group.smx(), ccs, n_refs):
    accept = ''
    if params.symop_threshold:
      if cc > params.symop_threshold:
        true_symops.append(smx)
        accept = '***'
    print '%10s %6d %.3f %s' % (smx, n_ref, cc, accept)

  if params.symop_threshold:
    from cctbx.sgtbx import space_group as sgtbx_space_group
    sg = sgtbx_space_group()
    for symop in true_symops:
      sg = sg.expand_smx(symop)
    for ltr in space_group.ltr():
      sg = sg.expand_ltr(ltr)
    sg_symbols = sg.match_tabulated_settings()
    print ''
    print 'Derived point group from symmetry operations: %s' % \
      sg_symbols.hermann_mauguin()
    print ''

  return

def offset_miller_indices(miller_indices, offset):
  from dials.array_family import flex
  return flex.miller_index(
    *[mi.iround() for mi in (miller_indices.as_vec3_double() + offset).parts()])

def test_P1_crystal_indexing(reflections, experiment, params):
  if params.grid_search_scope == 0:
    return

  from copy import deepcopy
  from dials.array_family import flex

  original_miller_indices = reflections['miller_index']

  crystal = experiment.crystal
  space_group = crystal.get_space_group()
  unit_cell = crystal.get_unit_cell()
  from cctbx.crystal import symmetry as crystal_symmetry
  cs = crystal_symmetry(unit_cell, space_group.type().lookup_symbol())

  from cctbx.miller import set as miller_set

  data = reflections['intensity.sum.value'] / \
         flex.sqrt(reflections['intensity.sum.variance'])

  ms = miller_set(cs, original_miller_indices)
  ms = ms.array(data)

  if params.d_min or params.d_max:
    ms = ms.resolution_filter(d_min=params.d_min, d_max=params.d_max)

  if params.asu:
    ms = ms.map_to_asu()

  print 'Checking HKL origin:'
  print ''
  print 'dH dK dL %6s %5s' % ('Nref', 'CC')

  g = params.grid_search_scope

  for h in range(-g, g + 1):
    for k in range(-g, g + 1):
      for l in range(-g, g + 1):

        for smx in ['-x,-y,-z']:
          reindexed = deepcopy(reflections)
          # hkl offset doubled as equivalent of h0 + 1, hI - 1
          miller_indices = offset_miller_indices(reflections['miller_index'],
                                                 (2 * h, 2 * k, 2* l))
          reindexed_miller_indices = sgtbx.change_of_basis_op(smx).apply(
            miller_indices)
          rms = miller_set(cs, reindexed_miller_indices)
          rms = rms.array(data)
          if params.d_min or params.d_max:
            rms = rms.resolution_filter(d_min=params.d_min, d_max=params.d_max)

          if params.asu:
            rms = rms.map_to_asu()

          intensity, intensity_rdx = rms.common_sets(ms)
          cc = intensity.correlation(intensity_rdx).coefficient()

          if cc > params.symop_threshold or (h == k == l == 0):
            print '%2d %2d %2d %6d %.3f' % \
              (h, k, l, intensity.size(), cc)

  print ''

  return

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
  reflections = reflections[0]

  test_P1_crystal_indexing(reflections, experiment, params)
  test_crystal_pointgroup_symmetry(reflections, experiment, params)

if __name__ == '__main__':
  import sys
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run(sys.argv[1:])
