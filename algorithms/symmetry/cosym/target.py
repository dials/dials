from __future__ import absolute_import, division, print_function

import logging
logger = logging.getLogger(__name__)

from dials.util import log

debug_handle = log.debug_handle(logger)
info_handle = log.info_handle(logger)

from cctbx.array_family import flex
from cctbx import sgtbx
from cctbx import miller

class Target(object):

  def __init__(self, miller_arrays, weights=None, min_pairs=None,
               lattice_group=None, dimensions=None, verbose=False):

    self._miller_arrays = miller_arrays
    self.verbose = verbose
    if weights is not None:
      assert weights in ('count', 'standard_error')
    self._weights = weights
    self._min_pairs = min_pairs

    miller_array_all = None
    lattice_ids = None
    space_group = None
    lattice_id = -1

    for intensities in miller_arrays:
      assert intensities.is_unique_set_under_symmetry()
      lattice_id += 1
      if space_group is None:
        space_group = intensities.space_group()
      else:
        assert intensities.space_group() == space_group

      ids = intensities.customized_copy(
        data=flex.double(intensities.size(), lattice_id), sigmas=None)
      assert ids.size() == intensities.size()
      if miller_array_all is None:
        miller_array_all = intensities
        lattice_ids = ids
      else:
        miller_array_all = miller_array_all.customized_copy(
          indices=miller_array_all.indices().concatenate(intensities.indices()),
          data=miller_array_all.data().concatenate(intensities.data()),
          sigmas=None)
        lattice_ids = lattice_ids.customized_copy(
          indices=lattice_ids.indices().concatenate(ids.indices()),
          data=lattice_ids.data().concatenate(ids.data()))
      assert miller_array_all.size() == lattice_ids.size()

    data = miller_array_all.customized_copy(anomalous_flag=False)
    cb_op_to_primitive = data.change_of_basis_op_to_primitive_setting()
    data = data.change_basis(cb_op_to_primitive).map_to_asu()

    resort = True
    if resort:
      order = flex.sort_permutation(lattice_ids.data())
      sorted_lattice_id = flex.select(lattice_ids.data(), order)
      sorted_data = data.data().select( order)
      sorted_indices = data.indices().select( order)
      self._lattice_ids = sorted_lattice_id
      self._data = data.customized_copy(indices = sorted_indices, data=sorted_data)
    else:
      self._lattice_ids = self._lattice_ids.data() # type flex int
      self._data = data # type miller array with flex double data
    assert isinstance(self._data.indices(), type(flex.miller_index()))
    assert isinstance(self._data.data(), type(flex.double()))

    # construct a lookup for the separate lattices
    last_id = -1; self._lattices = flex.int()
    for n in xrange(len(self._lattice_ids)):
      if self._lattice_ids[n] != last_id:
        last_id = self._lattice_ids[n]
        self._lattices.append(n)

    self._sym_ops = set(['x,y,z'])
    if lattice_group is None:
      self._sym_ops.update(set(
        [op.operator.as_xyz() for op in self.generate_twin_operators()]))
    else:
      self._sym_ops.update(set(
        [op.as_xyz() for op in lattice_group.group()]))
    if dimensions is None:
      dimensions = max(2, len(self._sym_ops))
    self.set_dimensions(dimensions)

    import copy
    self._lattice_group = copy.deepcopy(self._data.space_group())
    for sym_op in self._sym_ops:
      self._lattice_group.expand_smx(sym_op)

    logger.debug(
      'Lattice group: %s (%i symops)' %(
        self._lattice_group.info().symbol_and_number(), len(self._lattice_group)))

    import time
    t0 = time.time()
    self.compute_rij_wij()
    t1 = time.time()
    import scitbx.math
    rij = self.rij_matrix.as_1d()
    non_zero_sel = rij != 0
    logger.debug('Computed Rij matrix in %.2f seconds' %(t1 - t0))
    logger.debug('%i (%.1f%%) non-zero elements of Rij matrix' %(
        non_zero_sel.count(True), 100*non_zero_sel.count(True)/non_zero_sel.size()))
    scitbx.math.basic_statistics(rij.select(non_zero_sel)).show(f=debug_handle)

    return

  def set_dimensions(self, dimensions):
    self.dim = dimensions
    logger.info('Using %i dimensions for analysis' %self.dim)

  def generate_twin_operators(self):
    from mmtbx.scaling.twin_analyses import twin_laws
    TL = twin_laws(miller_array=self._data)
    for twin_law in TL.operators:
      if self.verbose: logger.info(twin_law.operator.r().as_hkl())
    return TL.operators

  def lattice_lower_upper_index(self, lattice_id):
    lower_index = self._lattices[lattice_id]
    upper_index = None
    if lattice_id < len(self._lattices)-1:
      upper_index = self._lattices[lattice_id+1]
    else:
      assert lattice_id == len(self._lattices)-1
    return lower_index, upper_index

  def compute_rij_wij(self, use_cache=True):

    group = flex.bool(self._lattices.size(), True)

    n_lattices = group.count(True)
    n_sym_ops = len(self._sym_ops)

    NN = n_lattices * n_sym_ops

    index_selected = group.iselection()
    self.rij_matrix = flex.double(flex.grid(NN,NN),0.)
    if self._weights is None:
      self.wij_matrix = None
    else:
      self.wij_matrix = flex.double(flex.grid(NN,NN),0.)

    indices = {}
    space_group_type = self._data.space_group().type()
    for cb_op in self._sym_ops:
      cb_op = sgtbx.change_of_basis_op(cb_op)
      indices_reindexed = cb_op.apply(self._data.indices())
      miller.map_to_asu(space_group_type, False, indices_reindexed)
      indices[cb_op.as_xyz()] = indices_reindexed

    rij_cache = {}

    for i in xrange(n_lattices):
      for j in xrange(i, n_lattices):
        for k, cb_op_k in enumerate(self._sym_ops):
          cb_op_k = sgtbx.change_of_basis_op(cb_op_k)
          for kk, cb_op_kk in enumerate(self._sym_ops):
            if i == j and k == kk:
              # don't include correlation of dataset with itself
              continue
            cb_op_kk = sgtbx.change_of_basis_op(cb_op_kk)

            ik = i + (n_lattices * k)
            jk = j + (n_lattices * kk)

            key = (i, j, str(cb_op_k.inverse() * cb_op_kk))
            if use_cache and key in rij_cache:
              cc, n = rij_cache[key]
            else:
              i_lower, i_upper = self.lattice_lower_upper_index(i)
              j_lower, j_upper = self.lattice_lower_upper_index(j)

              indices_i = indices[cb_op_k.as_xyz()][i_lower:i_upper]
              indices_j = indices[cb_op_kk.as_xyz()][j_lower:j_upper]
              intensities_i = self._data.data()[i_lower:i_upper]
              intensities_j = self._data.data()[j_lower:j_upper]

              matches = miller.match_indices(indices_i, indices_j)
              pairs = matches.pairs()
              isel_i = pairs.column(0)
              isel_j = pairs.column(1)
              isel_i = isel_i.select(self._lattice_group.epsilon(indices_i.select(isel_i)) == 1)
              isel_j = isel_j.select(self._lattice_group.epsilon(indices_j.select(isel_j)) == 1)
              corr = flex.linear_correlation(
                intensities_i.select(isel_i),
                intensities_j.select(isel_j))

              if corr.is_well_defined():
                cc = corr.coefficient()
                n = corr.n()
                rij_cache[key] = (cc, n)
              else:
                cc = None
                n = None

            if n < self._min_pairs:
              continue

            if cc is not None and n is not None:
              if self._weights == 'count':
                self.wij_matrix[(ik, jk)] = n
                self.wij_matrix[(jk, ik)] = n
              elif self._weights == 'standard_error':
                assert n > 2
                # http://www.sjsu.edu/faculty/gerstman/StatPrimer/correlation.pdf
                import math
                se = math.sqrt((1-cc**2)/(n-2))
                wij = 1/se
                self.wij_matrix[(ik, jk)] = wij
                self.wij_matrix[(jk, ik)] = wij

              self.rij_matrix[(ik, jk)] = cc
              self.rij_matrix[(jk, ik)] = cc

    return self.rij_matrix, self.wij_matrix

  def compute_functional(self, x):
    # x is a flattened list of the N-dimensional vectors, i.e. coordinates in
    # the first dimension are stored first, followed by the coordinates in the
    # second dimension, etc.

    assert (x.size() // self.dim) == (self._lattices.size() * len(self._sym_ops))
    inner = self.rij_matrix.deep_copy()
    NN = x.size() // self.dim
    for i in range(self.dim):
      coord = x[i*NN:(i+1)*NN]
      outer_prod = coord.matrix_outer_product(coord)
      inner -= outer_prod
    elements = inner*inner
    if self.wij_matrix is not None:
      elements = self.wij_matrix*elements
    f = 0.5 * flex.sum(elements)
    return f

  def compute_gradients_fd(self, x, eps=1e-6):
    grad = flex.double(x.size(), 0)
    for i in range(grad.size()):
      x[i] += eps # x + eps
      fp = self.compute_functional(x)
      x[i] -= (2 * eps) # x - eps
      fm = self.compute_functional(x)
      x[i] += eps # reset to original values
      grad[i] += ((fp - fm)/(2*eps))
    return grad

  def compute_functional_and_gradients(self, x):
    f = self.compute_functional(x)
    grad = flex.double()
    if self.wij_matrix is not None:
      wrij_matrix = self.wij_matrix * self.rij_matrix
    else:
      wrij_matrix = self.rij_matrix

    coords = []
    NN = x.size() // self.dim
    for i in range(self.dim):
      coords.append(x[i*NN:(i+1)*NN])

    # term 1
    for i in range(self.dim):
      grad.extend(wrij_matrix.matrix_multiply(coords[i]))

    for i in range(self.dim):
      tmp_array = flex.double()
      tmp = coords[i].matrix_outer_product(coords[i])
      if self.wij_matrix is not None:
        tmp = self.wij_matrix * tmp
      for j in range(self.dim):
        tmp_array.extend(tmp.matrix_multiply(coords[j]))
      grad -= tmp_array
    grad *= -2

    #grad_fd = self.compute_gradients_fd(x)
    #assert grad.all_approx_equal_relatively(grad_fd, relative_error=1e-4)

    return f, grad

  def curvatures(self, x):

    coords = []
    NN = x.size() // self.dim
    for i in range(self.dim):
      coords.append(x[i*NN:(i+1)*NN])

    curvs = flex.double()
    if self.wij_matrix is not None:
      wij = self.wij_matrix
    else:
      wij = flex.double(self.rij_matrix.accessor(), 1)
    for i in range(self.dim):
      curvs.extend(wij.matrix_multiply(coords[i] * coords[i]))
    curvs *= 2

    #curvs_fd = self.curvatures_fd(x)
    #assert curvs.all_approx_equal_relatively(curvs_fd, relative_error=1e-2)

    return curvs

  def curvatures_fd(self, x, eps=1e-6):
    f = self.compute_functional(x)
    curvs = flex.double(x.size(), 0)
    for i in range(curvs.size()):
      x[i] += eps # x + eps
      fp = self.compute_functional(x)
      x[i] -= (2 * eps) # x - eps
      fm = self.compute_functional(x)
      x[i] += eps # reset to original values
      curvs[i] += ((fm - 2 * f + fp)/(eps**2))
    return curvs

  def plot_rij_matrix(self, plot_name=None):
    if self.rij_matrix.all()[0] > 2000:
      return
    from matplotlib import pyplot as plt
    fig = plt.figure(figsize=(10,8))
    plt.clf()
    plt.imshow(self.rij_matrix.as_numpy_array(), vmin=-1, vmax=1, cmap='PiYG')
    plt.colorbar()
    if plot_name is not None:
      plt.savefig(plot_name)
    else:
      plt.show()

  def plot_wij_matrix(self, plot_name=None):
    if self._weights is None or self.wij_matrix.all()[0] > 2000:
      return
    from matplotlib import pyplot as plt
    fig = plt.figure(figsize=(10,8))
    plt.clf()
    plt.imshow(self.wij_matrix.as_numpy_array(), vmin=0, vmax=flex.max(self.wij_matrix))
    plt.colorbar()
    if plot_name is not None:
      plt.savefig(plot_name)
    else:
      plt.show()

  def plot_rij_histogram(self, plot_name=None):
    rij = self.rij_matrix.as_1d()
    rij = rij.select(rij != 0)
    hist = flex.histogram(rij, data_min=-1, data_max=1, n_slots=100)
    logger.debug('Histogram of Rij values:')
    hist.show(f=debug_handle)
    from matplotlib import pyplot as plt
    fig = plt.figure(figsize=(10,8))
    plt.clf()
    plt.bar(hist.slot_centers(), hist.slots(), width=hist.slot_width())
    fontsize = 24
    plt.xlabel(r'$r_{ij}$', size=fontsize)
    plt.ylabel('Frequency', size=fontsize)
    plt.tick_params(axis='both', which='both', labelsize=fontsize)
    plt.tight_layout()
    if plot_name is not None:
      plt.savefig(plot_name, dpi=300)
    else:
      plt.show()

  def plot_wij_histogram(self, plot_name=None):
    if self._weights is None:
      return
    wij = self.wij_matrix.as_1d()
    hist = flex.histogram(wij, n_slots=50)
    logger.debug('Histogram of Wij values:')
    hist.show(f=debug_handle)
    from matplotlib import pyplot as plt
    fig = plt.figure(figsize=(10,8))
    plt.clf()
    plt.bar(hist.slot_centers(), hist.slots(), width=hist.slot_width())
    plt.yscale('log')
    plt.xlabel(r'$w_{ij}$')
    plt.ylabel('Frequency')
    if plot_name is not None:
      plt.savefig(plot_name)
    else:
      plt.show()

  def plot_rij_cumulative_frequency(self, plot_name=None):
    rij = self.rij_matrix.as_1d()
    perm = flex.sort_permutation(rij)
    from matplotlib import pyplot as plt
    fig = plt.figure(figsize=(10,8))
    plt.clf()
    plt.plot(rij.select(perm), flex.int_range(perm.size()))
    plt.xlabel(r'$r_{ij}$')
    plt.ylabel('Cumulative requency')
    if plot_name is not None:
      plt.savefig(plot_name)
    else:
      plt.show()

  def plot_wij_cumulative_frequency(self, plot_name=None):
    if self._weights is None:
      return
    wij = self.wij_matrix.as_1d()
    perm = flex.sort_permutation(wij)
    import scitbx.math
    non_zero_sel = wij > 0
    logger.info('%i (%.1f%%) non-zero elements of Wij matrix' %(
      non_zero_sel.count(True), 100*non_zero_sel.count(True)/non_zero_sel.size()))
    scitbx.math.basic_statistics(wij.select(non_zero_sel)).show(f=debug_handle)
    from matplotlib import pyplot as plt
    fig = plt.figure(figsize=(10,8))
    plt.clf()
    plt.plot(wij.select(perm), flex.int_range(perm.size()))
    plt.xlabel(r'$w_{ij}$')
    plt.ylabel('Cumulative requency')
    if plot_name is not None:
      plt.savefig(plot_name)
    else:
      plt.show()

  def get_sym_ops(self):
    return self._sym_ops
