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

from __future__ import division

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
    refl['intensity.sum.value'] / flex.sqrt(refl['intensity.sum.variance']))

def offset_miller_indices(indices, offset):
  from dials.array_family import flex
  return flex.miller_index(
    *[mi.iround() for mi in (indices.as_vec3_double() + offset).parts()])

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
    if map_to_asu:
      rms = rms.map_to_asu()
      _ms = _ms.map_to_asu()
    intensity, intensity_rdx = rms.common_sets(_ms)
    cc = intensity.correlation(intensity_rdx).coefficient()
    ccs.append(cc)
    offsets.append(hkl)
    nref.append(intensity.size())

  return offsets, ccs, nref
