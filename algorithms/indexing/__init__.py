from __future__ import division

from cctbx.array_family import flex
from dials_algorithms_indexing_ext import *

def index_reflections(
    reflections, reciprocal_space_points, crystal_models, d_min,
    tolerance=0.3, verbosity=0):
  if 'miller_index' not in reflections:
    reflections['miller_index'] = flex.miller_index(len(reflections))
  if d_min is not None:
    d_spacings = 1/reciprocal_space_points.norms()
    inside_resolution_limit = d_spacings > d_min
  else:
    inside_resolution_limit = flex.bool(reciprocal_space_points.size(), True)
  sel = inside_resolution_limit & (reflections['id'] == -1)
  isel = sel.iselection()
  rlps = reciprocal_space_points.select(isel)
  refs = reflections.select(isel)

  diffs = []
  norms = []
  hkl_ints = []

  UB_matrices = flex.mat3_double([cm.get_A() for cm in crystal_models])

  if 1:
    # Use fast c++ version
    result = AssignIndices(rlps, UB_matrices, tolerance=tolerance)
    miller_indices = result.miller_indices()
    crystal_ids = result.crystal_ids()
    n_rejects = result.n_rejects()

    reflections['miller_index'].set_selected(sel, miller_indices)
    reflections['id'].set_selected(sel, crystal_ids)

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
    for i_lattice, crystal_model in enumerate(crystal_models):
      print "model %i (%i reflections):" %(
        i_lattice+1, (reflections['id'] == i_lattice).count(True))
      print crystal_model
      print

    print "%i unindexed reflections" %n_rejects


def index_reflections_local(
    reflections, reciprocal_space_points, crystal_models, d_min,
    tolerance=0.3, verbosity=0):
  if 'miller_index' not in reflections:
    reflections['miller_index'] = flex.miller_index(len(reflections))
  if d_min is not None:
    d_spacings = 1/reciprocal_space_points.norms()
    inside_resolution_limit = d_spacings > d_min
  else:
    inside_resolution_limit = flex.bool(reciprocal_space_points.size(), True)
  sel = inside_resolution_limit & (reflections['id'] == -1)
  isel = sel.iselection()
  rlps = reciprocal_space_points.select(isel)
  refs = reflections.select(isel)

  assert len(crystal_models) == 1

  diffs = []
  norms = []
  hkl_ints = []

  difference_vectors = flex.vec3_double()

  k = 10
  from annlib_ext import AnnAdaptor as ann_adaptor
  ann = ann_adaptor(data=rlps.as_double(), dim=3, k=k)
  ann.query(rlps.as_double())
  from scitbx import matrix
  #for i in range(rlps.size()):
    #r_i = matrix.col(rlps[i])
    #for ik in range(k):
      #r_ik = matrix.col(rlps[ann.nn[(i*k + ik)]])
      #d = r_i - r_ik
      #difference_vectors.append(d.elems)

  import math
  epsilon = 0.05
  delta = 5

  from libtbx.math_utils import nearest_integer as nint
  from libtbx.test_utils import approx_equal

  for i_lattice, crystal_model in enumerate(crystal_models):
    A = crystal_model.get_A()
    A_inv = A.inverse()

    import networkx as nx
    G = nx.Graph()
    G.add_nodes_from(range(len(rlps)))

    for i in range(rlps.size()):
      i_k = i * k
      for i_ann in range(k):
        i_k_plus_i_ann = i_k + i_ann
        j = ann.nn[i_k_plus_i_ann]
        if i > j and G.has_edge(i, j):
          continue
        d_r = matrix.col(rlps[i]) - matrix.col(rlps[j])
        h_f = A_inv * d_r
        h_ij = matrix.col([nint(f) for f in h_f.elems])
        d_h = h_f - h_ij
        n_h = d_h.length()
        l_ij = 1 - math.exp(-2 * sum(
          [(max(abs(d_h[ii]) - epsilon, 0)/epsilon)**2 +
           (max(abs(h_ij[ii]) - delta, 0))**2 for ii in range(3)]))

        if i < j:
          G.add_edge(i, j, weight=l_ij, h_ij=h_ij, d_r=d_r, h_f=h_f)
        else:
          G.add_edge(j, i, weight=l_ij, h_ij=-h_ij, d_r=-d_r, h_f=-h_f)


    T = nx.minimum_spanning_tree(G)

    hkl = flex.vec3_int(T.number_of_nodes())
    subtree_id = flex.int(T.number_of_nodes(), -1)
    last_hkl = (0,0,0)
    subtree_id[0] = 0
    l_min = 0.8 # XDS parameter INDEX_QUALITY=

    for i, j in nx.dfs_edges(T, source=0):
      if j < i:
        sign = -1
      else:
        sign = 1
      hkl_i = hkl[i]
      edge = T.get_edge_data(i, j)
      l_ij = edge['weight']
      h_ij = sign * edge['h_ij']
      hkl[j] = (matrix.col(hkl_i) - matrix.col(h_ij)).elems
      if l_ij < l_min:
        subtree_id[j] = subtree_id[i]
      else:
        subtree_id[j] = subtree_id[i] + 1

    unique_subtree_ids = set(subtree_id)
    largest_subtree_id = -1
    largest_subtree_size = 0
    for i in unique_subtree_ids:
      subtree_size = subtree_id.count(i)
      if subtree_size > largest_subtree_size:
        largest_subtree_size = subtree_size
        largest_subtree_id = i

    print "Largest subtree: %i, %i" %(largest_subtree_id, largest_subtree_size)

    subtree_sel = (subtree_id == largest_subtree_id)

  n_rejects = 0

  ref_sel = refs.select(subtree_sel)
  rlp_sel = rlps.select(subtree_sel)
  hkl_sel = hkl.select(subtree_sel).as_vec3_double()

  # need to do a proper search to find best offset
  hf_0 = A_inv * rlp_sel[0]
  h_0 = matrix.col([nint(i) for i in hf_0.elems])
  offset = (h_0 - matrix.col(hkl_sel[0])).elems

  h = hkl_sel + flex.vec3_double(hkl_sel.size(), offset)
  test_rlp = tuple(A) * h
  d_rlp = rlp_sel - test_rlp
  print d_rlp.min(), d_rlp.max(), d_rlp.mean()
  for i in range(20):
    print "(%i, %i, %i)" %h[i], "(%.2f, %.2f, %.2f)" %test_rlp[i], "(%.2f, %.2f, %.2f)" %d_rlp[i]

  refs['miller_index'].set_selected(
    subtree_sel, flex.miller_index(list(h.iround())))
  refs['id'].set_selected(subtree_sel, 0)

  reflections['miller_index'].set_selected(isel, refs['miller_index'])
  reflections['id'].set_selected(isel, refs['id'])

  if verbosity > 0:
    for i_lattice, crystal_model in enumerate(crystal_models):
      print "model %i (%i reflections):" %(
        i_lattice+1, (reflections['id'] == i_lattice).count(True))
      print crystal_model
      print

    print "%i unindexed reflections" %n_rejects

#index_reflections = index_reflections_local
