from __future__ import division

from cctbx.array_family import flex

def index_reflections(
    reflections, reciprocal_space_points, crystal_models, d_min,
    tolerance=0.3, verbosity=0):

  if d_min is not None:
    d_spacings = 1/reciprocal_space_points.norms()
    inside_resolution_limit = d_spacings > d_min
  else:
    inside_resolution_limit = flex.bool(reciprocal_space_points.size(), True)
  sel = inside_resolution_limit & (reflections.crystal() == -1)
  isel = sel.iselection()
  rlps = reciprocal_space_points.select(isel)

  diffs = []
  norms = []
  hkl_ints = []

  for i_lattice, crystal_model in enumerate(crystal_models):
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
    reflections[i_ref].miller_index = miller_index
    reflections[i_ref].crystal = i_best_lattice

  if verbosity > 0:
    for i_lattice, crystal_model in enumerate(crystal_models):
      print "model %i (%i reflections):" %(
        i_lattice+1, (reflections.crystal() == i_lattice).count(True))
      print crystal_model
      print

    print "%i unindexed reflections" %n_rejects
