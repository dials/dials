from __future__ import absolute_import, division, print_function

import logging
logger = logging.getLogger(__name__)

import copy
import math

from dials.util import log

debug_handle = log.debug_handle(logger)
info_handle = log.info_handle(logger)

from cctbx.array_family import flex
from cctbx import sgtbx
from cctbx import miller
import cctbx.sgtbx.cosets

class Target(object):
  """Target function for cosym analysis.

  Attributes:
    dim (int): The number of dimensions used in the analysis.

  """

  def __init__(self, intensities, lattice_ids, weights=None, min_pairs=None,
               lattice_group=None, dimensions=None, nproc=1):
    r"""Intialise a Target object.

    Args:
      intensities (cctbx.miller.array): The intensities on which to perform
        cosym anaylsis.
      lattice_ids (scitbx.array_family.flex.int): An array of equal size to
        `intensities` which maps each reflection to a given lattice (dataset).
      weights (str): Optionally include weights in the target function.
        Allowed values are `None`, "count" and "standard_error". The default
        is to use no weights. If "count" is set, then weights are equal to the
        number of pairs of reflections used in calculating each value of the
        rij matrix. If "standard_error" is used, then weights are defined as
        :math:`w_{ij} = 1/s`, where :math:`s = \sqrt{(1-r_{ij}^2)/(n-2)}`.
        See also http://www.sjsu.edu/faculty/gerstman/StatPrimer/correlation.pdf.
      min_pairs (int): Only calculate the correlation coefficient between two
        datasets if they have more than `min_pairs` of common reflections.
      lattice_group (cctbx.sgtbx.space_group): Optionally set the lattice
        group to be used in the analysis.
      dimensions (int): Optionally override the number of dimensions to be used
        in the analysis. If not set, then the number of dimensions used is
        equal to the greater of 2 or the number of symmetry operations in the
        lattice group.
      nproc (int): number of processors to use for computing the rij matrix.

    """
    if weights is not None:
      assert weights in ('count', 'standard_error')
    self._weights = weights
    self._min_pairs = min_pairs
    self._nproc = nproc

    data = intensities.customized_copy(anomalous_flag=False)
    cb_op_to_primitive = data.change_of_basis_op_to_primitive_setting()
    data = data.change_basis(cb_op_to_primitive).map_to_asu()

    order = flex.sort_permutation(lattice_ids)
    sorted_lattice_id = flex.select(lattice_ids, order)
    sorted_data = data.data().select( order)
    sorted_indices = data.indices().select( order)
    self._lattice_ids = sorted_lattice_id
    self._data = data.customized_copy(indices=sorted_indices, data=sorted_data)
    assert isinstance(self._data.indices(), type(flex.miller_index()))
    assert isinstance(self._data.data(), type(flex.double()))

    # construct a lookup for the separate lattices
    last_id = -1; self._lattices = flex.int()
    for n in xrange(len(self._lattice_ids)):
      if self._lattice_ids[n] != last_id:
        last_id = self._lattice_ids[n]
        self._lattices.append(n)

    self._sym_ops = set(['x,y,z'])
    self._lattice_group = lattice_group
    self._sym_ops.update(
      set([op.as_xyz() for op in self._generate_twin_operators()]))
    if dimensions is None:
      dimensions = max(2, len(self._sym_ops))
    self.set_dimensions(dimensions)

    self._lattice_group = copy.deepcopy(self._data.space_group())
    for sym_op in self._sym_ops:
      self._lattice_group.expand_smx(sym_op)
    self._patterson_group = self._lattice_group.build_derived_patterson_group()

    logger.debug(
      'Lattice group: %s (%i symops)' %(
        self._lattice_group.info().symbol_and_number(), len(self._lattice_group)))
    logger.debug(
      'Patterson group: %s' %self._patterson_group.info().symbol_and_number())

    self._compute_rij_wij()

  def set_dimensions(self, dimensions):
    """Set the number of dimensions for analysis.

    Args:
      dimensions (int): The number of dimensions to be used.

    """
    self.dim = dimensions
    logger.info('Using %i dimensions for analysis' %self.dim)

  def _generate_twin_operators(self, lattice_symmetry_max_delta=3.0):
    # see also mmtbx.scaling.twin_analyses.twin_laws
    cb_op_to_niggli_cell = \
      self._data.change_of_basis_op_to_niggli_cell()
    if self._lattice_group is None:
      minimum_cell_symmetry = self._data.crystal_symmetry().change_basis(
        cb_op=cb_op_to_niggli_cell)
      self._lattice_group = sgtbx.lattice_symmetry.group(
        reduced_cell=minimum_cell_symmetry.unit_cell(),
        max_delta=lattice_symmetry_max_delta)
      intensity_symmetry = minimum_cell_symmetry.reflection_intensity_symmetry(
        anomalous_flag=self._data.anomalous_flag())
      cb_op = cb_op_to_niggli_cell.inverse()
    else:
      cb_op = sgtbx.change_of_basis_op()
      intensity_symmetry = self._data.reflection_intensity_symmetry()

    operators = []
    for partition in sgtbx.cosets.left_decomposition(
      g = self._lattice_group,
      h = intensity_symmetry.space_group()
          .build_derived_acentric_group()
          .make_tidy()).partitions[1:] :
      if (partition[0].r().determinant() > 0):
        operators.append(cb_op.apply(partition[0]))

    return operators

  def _lattice_lower_upper_index(self, lattice_id):
    lower_index = self._lattices[lattice_id]
    upper_index = None
    if lattice_id < len(self._lattices)-1:
      upper_index = self._lattices[lattice_id+1]
    else:
      assert lattice_id == len(self._lattices)-1
    return lower_index, upper_index

  def _compute_rij_wij(self, use_cache=True):
    """Compute the rij_wij matrix."""
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

    def _compute_rij_matrix_one_row_block(i):
      rij_cache = {}

      n_sym_ops = len(self._sym_ops)
      NN = n_lattices * n_sym_ops

      from scipy import sparse
      rij_row = []
      rij_col = []
      rij_data = []
      if self._weights is not None:
        wij_row = []
        wij_col = []
        wij_data = []
      else:
        wij = None

      i_lower, i_upper = self._lattice_lower_upper_index(i)
      intensities_i = self._data.data()[i_lower:i_upper]

      for j in range(i, n_lattices):

        j_lower, j_upper = self._lattice_lower_upper_index(j)
        intensities_j = self._data.data()[j_lower:j_upper]

        for k, cb_op_k in enumerate(self._sym_ops):
          cb_op_k = sgtbx.change_of_basis_op(cb_op_k)

          indices_i = indices[cb_op_k.as_xyz()][i_lower:i_upper]

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
              indices_j = indices[cb_op_kk.as_xyz()][j_lower:j_upper]

              matches = miller.match_indices(indices_i, indices_j)
              pairs = matches.pairs()
              isel_i = pairs.column(0)
              isel_j = pairs.column(1)
              isel_i = isel_i.select(
                self._patterson_group.epsilon(indices_i.select(isel_i)) == 1)
              isel_j = isel_j.select(
                self._patterson_group.epsilon(indices_j.select(isel_j)) == 1)
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
                wij_row.extend([ik, jk])
                wij_col.extend([jk, ik])
                wij_data.extend([n, n])
              elif self._weights == 'standard_error':
                assert n > 2
                # http://www.sjsu.edu/faculty/gerstman/StatPrimer/correlation.pdf
                se = math.sqrt((1-cc**2)/(n-2))
                wij = 1/se
                wij_row.extend([ik, jk])
                wij_col.extend([jk, ik])
                wij_data.extend([wij, wij])

              rij_row.extend([ik, jk])
              rij_col.extend([jk, ik])
              rij_data.extend([cc, cc])

      rij = sparse.coo_matrix((rij_data, (rij_row, rij_col)), shape=(NN, NN))
      if self._weights is not None:
        wij = sparse.coo_matrix((wij_data, (wij_row, wij_col)), shape=(NN, NN))

      return rij, wij

    from libtbx import easy_mp
    args = [(i,) for i in range(n_lattices)]
    results = easy_mp.parallel_map(
      _compute_rij_matrix_one_row_block,
      args,
      processes=self._nproc,
      iterable_type=easy_mp.posiargs,
      method='multiprocessing')

    rij_matrix = None
    wij_matrix = None
    for i, (rij, wij) in enumerate(results):
      if rij_matrix is None:
        rij_matrix = rij
      else:
        rij_matrix += rij
      if wij is not None:
        if wij_matrix is None:
          wij_matrix = wij
        else:
          wij_matrix += wij

    self.rij_matrix = flex.double(rij_matrix.todense())
    if wij_matrix is not None:
      import numpy as np
      self.wij_matrix = flex.double(wij_matrix.todense().astype(np.float64))

    return self.rij_matrix, self.wij_matrix

  def compute_functional(self, x):
    """Compute the target function at coordinates `x`.

    Args:
      x (scitbx.array_family.flex.double):
        a flattened list of the N-dimensional vectors, i.e. coordinates in
        the first dimension are stored first, followed by the coordinates in
        the second dimension, etc.

    Returns:
      f (float): The value of the target function at coordinates `x`.

    """
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
    """Compute the gradients at coordinates `x` using finite differences.

    Args:
      x (scitbx.array_family.flex.double):
        a flattened list of the N-dimensional vectors, i.e. coordinates in
        the first dimension are stored first, followed by the coordinates in
        the second dimension, etc.
      eps (float):
        The value of epsilon to use in finite difference calculations.

    Returns:
      grad (scitbx.array_family.flex.double):
      The gradients of the target function with respect to the parameters.

    """
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
    """Compute the target function and gradients at coordinates `x`.

    Args:
      x (scitbx.array_family.flex.double):
        a flattened list of the N-dimensional vectors, i.e. coordinates in
        the first dimension are stored first, followed by the coordinates in
        the second dimension, etc.

    Returns:
      Tuple[float, scitbx.array_family.flex.double]:
      f: The value of the target function at coordinates `x`.
      grad: The gradients of the target function with respect to the parameters.

    """
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
    """Compute the curvature of the target function.

    Args:
      x (scitbx.array_family.flex.double):
        a flattened list of the N-dimensional vectors, i.e. coordinates in
        the first dimension are stored first, followed by the coordinates in
        the second dimension, etc.

    Returns:
      curvs (scitbx.array_family.flex.double):
      The curvature of the target function with respect to the parameters.

    """
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

    return curvs

  def curvatures_fd(self, x, eps=1e-6):
    """Compute the curvatures at coordinates `x` using finite differences.

    Args:
      x (scitbx.array_family.flex.double):
        a flattened list of the N-dimensional vectors, i.e. coordinates in
        the first dimension are stored first, followed by the coordinates in
        the second dimension, etc.
      eps (float):
        The value of epsilon to use in finite difference calculations.

    Returns:
      curvs (scitbx.array_family.flex.double):
      The curvature of the target function with respect to the parameters.

    """
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
    """Plot the matrix of rij values.

    Args:
      plot_name (str): The file name to save the plot to.
        If this is not defined then the plot is displayed in interactive mode.

    """
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
    """Plot the matrix of wij values.

    Args:
      plot_name (str): The file name to save the plot to.
        If this is not defined then the plot is displayed in interactive mode.

    """
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
    """Plot a histogram of the rij values.

    Args:
      plot_name (str): The file name to save the plot to.
        If this is not defined then the plot is displayed in interactive mode.

    """
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
    """Plot a histogram of the wij values.

    Args:
      plot_name (str): The file name to save the plot to.
        If this is not defined then the plot is displayed in interactive mode.

    """
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
    """Plot the cumulative frequency of the rij values.

    Args:
      plot_name (str): The file name to save the plot to.
        If this is not defined then the plot is displayed in interactive mode.

    """
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
    """Plot the cumulative frequency of the wij values.

    Args:
      plot_name (str): The file name to save the plot to.
        If this is not defined then the plot is displayed in interactive mode.

    """
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
    """Get the list of symmetry operations used in the analysis.

    Returns:
      List[cctbx.sgtbx.rt_mx]: The list of symmetry operations.

    """
    return self._sym_ops
