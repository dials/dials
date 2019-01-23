#!/usr/bin/env dials.python
#
# dials.algorithms.symmetry.origin
#
#  Copyright (C) 2016 Diamond Light Source
#
#  Author: Graeme Winter
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#
# Analysis of the origin of the diffraction pattern based on indexed and
# measured intensities.

from __future__ import absolute_import, division

def cctbx_crystal_from_dials(crystal):
  space_group = crystal.get_space_group()
  unit_cell = crystal.get_unit_cell()
  from cctbx.crystal import symmetry as crystal_symmetry
  return crystal_symmetry(unit_cell, space_group.type().lookup_symbol())

def cctbx_i_over_sigi_ms_from_dials_data(reflections, cctbx_crystal_symmetry):
  from dials.array_family import flex
  from cctbx.miller import set as miller_set
  refl = reflections.select(reflections['intensity.sum.variance'] > 0)
  return miller_set(cctbx_crystal_symmetry, refl['miller_index']).array(
    data=refl['intensity.sum.value'],
    sigmas=flex.sqrt(refl['intensity.sum.variance']))

def offset_miller_indices(indices, offset):
  from dials.array_family import flex
  return flex.miller_index(
    *[mi.iround() for mi in (indices.as_vec3_double() + offset).parts()])

def compute_miller_set_correlation(ms_a, ms_b, map_to_asu=False,
                                   merge_equivalents=False):
  """Compute correlation between two miller arrays.

  Args:
    ms_a (cctbx.miller.array): Input miller.array `a`.
    ms_b (cctbx.miller.array): Input miller.array `b`.
    map_to_asu (bool): If ``True``, then map miller indices to the asymmetric
      unit before matching miller indices between input miller arrays.
    merge_equivalents (bool): If ``True`` then merge symmetry equivalent
      reflections before matching miller indices between input miller arrays.

  Returns:
    tuple[int, float]: A tuple of the number of observations and the correlation
    coefficient.

  """
  if map_to_asu:
    # not obvious that this will help for the reasons stated below
    ms_a = ms_a.map_to_asu()
    ms_b = ms_b.map_to_asu()

  if merge_equivalents:
    # only want to do this if we have essentially "scaled" the data - if not
    # then we will get a smooth Wilson plot and about CC=1 (due to general
    # fall off with resolution)
    ms_a = ms_a.merge_equivalents().array()
    ms_b = ms_b.merge_equivalents().array()

  common_a, common_b = ms_a.common_sets(ms_b)

  return common_a.size(), common_a.correlation(common_b).coefficient()

def get_hkl_offset_correlation_coefficients(
  dials_reflections, dials_crystal, map_to_asu=False,
  grid_h=0, grid_k=0, grid_l=0, reference=None):

  # N.B. deliberately ignoring d_min, d_max as these are inconsistent with
  # changing the miller indices

  from dials.array_family import flex
  from cctbx.miller import set as miller_set
  from cctbx import sgtbx

  cs = cctbx_crystal_from_dials(dials_crystal)
  ms = cctbx_i_over_sigi_ms_from_dials_data(dials_reflections, cs)

  if reference:
    reference_ms = cctbx_i_over_sigi_ms_from_dials_data(reference, cs)
  else:
    reference_ms = None

  ccs = flex.double()
  offsets = flex.vec3_int()
  nref = flex.size_t()

  if reference:
    cb_op = sgtbx.change_of_basis_op('x,y,z')
  else:
    cb_op = sgtbx.change_of_basis_op('-x,-y,-z')

  hkl_test = [(h, k, l) for h in range(-grid_h, grid_h + 1) \
                        for k in range(-grid_k, grid_k + 1) \
                        for l in range(-grid_l, grid_l + 1)]

  for hkl in hkl_test:
    indices = offset_miller_indices(ms.indices(), hkl)
    reindexed_indices = cb_op.apply(indices)
    rms = miller_set(cs, reindexed_indices).array(ms.data())
    if reference_ms:
      _ms = reference_ms
    else:
      _ms = miller_set(cs, indices).array(ms.data())
    n, cc = compute_miller_set_correlation(_ms, rms, map_to_asu=map_to_asu)
    ccs.append(cc)
    offsets.append(hkl)
    nref.append(n)

  return offsets, ccs, nref
