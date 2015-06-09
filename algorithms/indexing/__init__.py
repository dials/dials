from __future__ import division

from cctbx.array_family import flex
from dials_algorithms_indexing_ext import *
from logging import info, debug

def index_reflections(
    reflections, experiments, d_min=None,
    tolerance=0.3, verbosity=0):
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

  if 1:
    # Use fast c++ version
    imgset_ids = reflections['imageset_id'].select(sel)

    for i_imgset, imgset in enumerate(experiments.imagesets()):
      sel_imgset = (imgset_ids == i_imgset)

      result = AssignIndices(
        rlps.select(sel_imgset), phi.select(sel_imgset), UB_matrices, tolerance=tolerance)

      miller_indices = result.miller_indices()
      crystal_ids = result.crystal_ids()
      n_rejects = result.n_rejects()

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

  else:
    # XXX Use old python version
    for i_lattice, crystal_model in enumerate(crystal_models):
      if crystal_model.num_scan_points > 0:
        from scitbx import matrix
        import math
        A_inv = [matrix.sqr(crystal.get_A_at_scan_point(
          math.floor(refs[i]['xyzobs.px.value'][2]))).inverse()
                 for i in range(len(refs))]
        hkl_float = flex.vec3_double(
          [A_inv[i] * rlps[i] for i in range(len(rlps))])
      else:
        A = crystal_model.get_A()
        A_inv = A.inverse()
        hkl_float = tuple(A_inv) * rlps
      hkl_int = hkl_float.iround()
      differences = hkl_float - hkl_int.as_vec3_double()

      diffs.append(differences)
      norms.append(differences.norms())
      hkl_ints.append(hkl_int)

    n_rejects = 0

    for i_hkl in range(hkl_int.size()):
      d = flex.vec3_double([diffs[j][i_hkl]
                       for j in range(len(crystal_models))])
      n = flex.double([norms[j][i_hkl]
                       for j in range(len(crystal_models))])
      potential_hkls = [hkl_ints[j][i_hkl]
                        for j in range(len(crystal_models))]

      i_best_lattice = flex.min_index(n)
      if n[i_best_lattice] > tolerance:
        n_rejects += 1
        continue
      miller_index = potential_hkls[i_best_lattice]
      i_ref = isel[i_hkl]
      reflections['miller_index'][i_ref] = miller_index
      reflections['id'][i_ref] = i_best_lattice

    # if more than one spot can be assigned the same miller index then choose
    # the closest one
    miller_indices = reflections['miller_index'].select(isel)
    for i_hkl, hkl in enumerate(miller_indices):
      if hkl == (0,0,0): continue
      iselection = (miller_indices == hkl).iselection()
      if len(iselection) > 1:
        for i in iselection:
          for j in iselection:
            if j <= i: continue
            crystal_i = reflections['id'][isel[i]]
            crystal_j = reflections['id'][isel[j]]
            if crystal_i != crystal_j:
              continue
            elif crystal_i == -1:
              continue
            assert hkl_ints[crystal_j][j] == hkl_ints[crystal_i][i]
            if norms[crystal_j][j] < norms[crystal_i][i]:
              i_ref = isel[i]
            else:
              i_ref = isel[j]
            reflections['miller_index'][i_ref] = (0,0,0)
            reflections['id'][i_ref] = -1

  if verbosity > 0:
    for i_expt, expt in enumerate(experiments):
      info("model %i (%i reflections):" %(
        i_expt+1, (reflections['id'] == i_expt).count(True)))
      info(expt.crystal)
      info("")

    info("%i unindexed reflections" %n_rejects)


def index_reflections_local(
    reflections, experiments, d_min=None,
    epsilon=0.05, delta=8, l_min=0.8, nearest_neighbours=20, verbosity=0):
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
  n_rejects = result.n_rejects()
  hkl = miller_indices.as_vec3_double().iround()

  n_rejects = (crystal_ids < 0).count(True)
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

  if verbosity > 0:
    for i_cryst, cryst in enumerate(experiments.crystals()):
      info("model %i (%i reflections):" %(
        i_cryst+1, (reflections['id'] == i_cryst).count(True)))
      info(cryst)
      info("")

    info("%i unindexed reflections" %n_rejects)
