from __future__ import division

from cctbx.array_family import flex
from dials_algorithms_indexing_ext import *

def index_reflections(reflections, experiments, d_min=None, tolerance=0.3):
  reciprocal_lattice_points = reflections['rlp']
  reflections['miller_index'] = flex.miller_index(len(reflections), (0,0,0))
  if d_min is not None:
    d_spacings = 1/reciprocal_lattice_points.norms()
    inside_resolution_limit = d_spacings > d_min
  else:
    inside_resolution_limit = flex.bool(reciprocal_lattice_points.size(), True)
  sel = inside_resolution_limit & (reflections['id'] == -1)
  isel = sel.iselection()
  rlps = reciprocal_lattice_points.select(isel)
  refs = reflections.select(isel)
  phi = refs['xyzobs.mm.value'].parts()[2]

  diffs = []
  norms = []
  hkl_ints = []

  UB_matrices = flex.mat3_double([cm.get_A() for cm in experiments.crystals()])
  imgset_ids = reflections['imageset_id'].select(sel)

  for i_imgset, imgset in enumerate(experiments.imagesets()):
    sel_imgset = (imgset_ids == i_imgset)

    result = AssignIndices(
      rlps.select(sel_imgset), phi.select(sel_imgset), UB_matrices, tolerance=tolerance)

    miller_indices = result.miller_indices()
    crystal_ids = result.crystal_ids()

    expt_ids = flex.int(crystal_ids.size(), -1)
    for i_cryst, cryst in enumerate(experiments.crystals()):
      sel_cryst = (crystal_ids == i_cryst)
      for i_expt in experiments.where(
        crystal=cryst, imageset=imgset):
        expt_ids.set_selected(sel_cryst, i_expt)

    reflections['miller_index'].set_selected(isel.select(sel_imgset), miller_indices)
    reflections['id'].set_selected(isel.select(sel_imgset), expt_ids)
    reflections.set_flags(
      reflections['miller_index'] != (0,0,0), reflections.flags.indexed)


def index_reflections_local(
    reflections, experiments, d_min=None,
    epsilon=0.05, delta=8, l_min=0.8, nearest_neighbours=20):
  from scitbx import matrix
  from libtbx.math_utils import nearest_integer as nint
  reciprocal_lattice_points = reflections['rlp']
  if 'miller_index' not in reflections:
    reflections['miller_index'] = flex.miller_index(len(reflections))
  if d_min is not None:
    d_spacings = 1/reciprocal_lattice_points.norms()
    inside_resolution_limit = d_spacings > d_min
  else:
    inside_resolution_limit = flex.bool(reciprocal_lattice_points.size(), True)
  sel = inside_resolution_limit & (reflections['id'] == -1)
  isel = sel.iselection()
  rlps = reciprocal_lattice_points.select(isel)
  refs = reflections.select(isel)
  phi = refs['xyzobs.mm.value'].parts()[2]

  diffs = []
  norms = []
  hkl_ints = []

  UB_matrices = flex.mat3_double([cm.get_A() for cm in experiments.crystals()])

  result = AssignIndicesLocal(
    rlps, phi, UB_matrices, epsilon=epsilon, delta=delta, l_min=l_min,
    nearest_neighbours=nearest_neighbours)
  miller_indices = result.miller_indices()
  crystal_ids = result.crystal_ids()
  hkl = miller_indices.as_vec3_double().iround()

  assert miller_indices.select(crystal_ids < 0).all_eq((0,0,0))

  for i_cryst in set(crystal_ids):
    if i_cryst < 0: continue

    A = experiments[i_cryst].crystal.get_A()
    A_inv = A.inverse()

    cryst_sel = crystal_ids == i_cryst
    ref_sel = refs.select(cryst_sel)
    rlp_sel = rlps.select(cryst_sel)
    hkl_sel = hkl.select(cryst_sel).as_vec3_double()

    d_sel = 1/rlp_sel.norms()
    d_perm = flex.sort_permutation(d_sel, reverse=True)

    hf_0 = A_inv * rlp_sel[d_perm[0]]
    h_0 = matrix.col([nint(j) for j in hf_0.elems])
    offset = h_0 - matrix.col(hkl_sel[d_perm[0]])
    #print "offset:", offset.elems

    h = hkl_sel + flex.vec3_double(hkl_sel.size(), offset.elems)

    refs['miller_index'].set_selected(
      cryst_sel, flex.miller_index(list(h.iround())))
    refs['id'].set_selected(cryst_sel, i_cryst)

  crystal_ids.set_selected(crystal_ids < 0, -1)
  refs['id'] = crystal_ids
  refs['miller_index'].set_selected(crystal_ids < 0, (0,0,0))

  reflections['miller_index'].set_selected(isel, refs['miller_index'])
  reflections['id'].set_selected(isel, refs['id'])
  reflections.set_flags(
    reflections['miller_index'] != (0,0,0), reflections.flags.indexed)
