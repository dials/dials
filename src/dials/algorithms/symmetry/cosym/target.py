"""Target function for cosym analysis."""

from __future__ import annotations

import concurrent.futures
import copy
import logging

import numpy as np
from ordered_set import OrderedSet
from scipy import sparse

import cctbx.sgtbx.cosets
from cctbx import miller, sgtbx
from cctbx.array_family import flex

from dials.algorithms.scaling.scaling_library import ExtendedDatasetStatistics
from dials_cosym_ext import matcher

logger = logging.getLogger(__name__)


def _lattice_lower_upper_index(lattices, lattice_id):
    lower_index = int(lattices[lattice_id])
    upper_index = None
    if lattice_id < len(lattices) - 1:
        upper_index = int(lattices[lattice_id + 1])
    else:
        assert lattice_id == len(lattices) - 1
    return lower_index, upper_index


class FakeArray:
    ## Confroms to a subset of the miller array interface.
    def __init__(self, data, sigmas):
        self._data = data
        self._sigmas = sigmas

    def data(self):
        return self._data

    def sigmas(self):
        return self._sigmas

    def size(self):
        return self._sigmas.size()


def _compute_rij_matrix_one_row_block(
    i,
    lattices,
    data,
    sym_ops,
    patterson_group,
    weights=True,  # use sigma weights or not for CC1/2 calculation
    min_pairs=3,
):
    ## Calculate the upper half-triangle of the rij matrix for row-block i.

    n_lattices = len(lattices)
    space_group_type = data.space_group().type()
    NN = n_lattices * len(sym_ops)

    ## Lists to aggregate data for matrix output.
    rij_row = []
    rij_col = []
    rij_data = []
    wij_row = []
    wij_col = []
    wij_data = []

    i_lower, i_upper = _lattice_lower_upper_index(lattices, i)
    intensities_i = data.data()[i_lower:i_upper]
    sigmas_i = data.sigmas()[i_lower:i_upper]
    cb_ops = [sgtbx.change_of_basis_op(cb_op_k) for cb_op_k in sym_ops]
    original_indices_i = data.indices()[i_lower:i_upper]

    for j in range(i, n_lattices):  # only calculating a half-triangle of rij
        j_lower, j_upper = _lattice_lower_upper_index(lattices, j)
        intensities_j = data.data()[j_lower:j_upper]
        sigmas_j = data.sigmas()[j_lower:j_upper]
        original_indices_j = data.indices()[j_lower:j_upper]
        rij_cache = {}

        for k, cb_op_k in enumerate(cb_ops):
            ## We initialise the miller index matcher per k, which creates a lookup map
            ## for dataset i reindexed by operator k. Then it is quicker to generate the
            ## matching indices in repeated lookups. Delay the construction until we
            ## require it, as if the data is in the cache we might not need it and it is
            ## relatively expensive.
            matcher_k = None
            cb_op_k_inverse = cb_op_k.inverse()
            for kk, cb_op_kk in enumerate(cb_ops):
                if i == j and k <= kk:
                    # don't include correlation of dataset with itself (i==j, k==kk)
                    # also make sure we're only filling a half-triangle for i==j
                    continue

                ik = i + (n_lattices * k)
                jk = j + (n_lattices * kk)

                key = str(cb_op_k_inverse * cb_op_kk)
                if key in rij_cache:
                    cc, n, n_pairs = rij_cache[key]
                else:
                    # only do expensive calculations at this point, as we may have
                    # been able to use the cache for the given combination of cb ops.
                    if not matcher_k:
                        indices_i = cb_op_k.apply(original_indices_i)
                        miller.map_to_asu(space_group_type, False, indices_i)
                        matcher_k = matcher(indices_i, patterson_group)

                    indices_j = cb_op_kk.apply(original_indices_j)
                    miller.map_to_asu(space_group_type, False, indices_j)
                    isel_i, isel_j = matcher_k.match(indices_j)

                    ma_j = FakeArray(
                        intensities_j.select(isel_j), sigmas_j.select(isel_j)
                    )
                    ma_i = FakeArray(
                        intensities_i.select(isel_i), sigmas_i.select(isel_i)
                    )
                    n_pairs = ma_i.size()
                    if ma_i.size() < min_pairs:
                        n, cc = (None, None)
                    else:
                        if weights:
                            corr, neff = ExtendedDatasetStatistics.weighted_cchalf(
                                ma_i, ma_j, assume_index_matching=True
                            )[0]
                            if neff:
                                cc = corr
                                n = neff
                            else:
                                n, cc = (None, None)
                        else:
                            cc = flex.linear_correlation(
                                ma_i.data(), ma_j.data()
                            ).coefficient()
                            n = n_pairs

                    rij_cache[key] = (cc, n, n_pairs)

                if (
                    n is None
                    or cc is None
                    or (min_pairs is not None and n_pairs < min_pairs)
                ):
                    continue

                wij_row.append(ik)
                wij_col.append(jk)
                wij_data.append(n)
                rij_row.append(ik)
                rij_col.append(jk)
                rij_data.append(cc)

    rij = sparse.coo_matrix((rij_data, (rij_row, rij_col)), shape=(NN, NN))
    wij = sparse.coo_matrix((wij_data, (wij_row, wij_col)), shape=(NN, NN))
    return rij, wij


class Target:
    """Target function for cosym analysis.

    Attributes:
      dim (int): The number of dimensions used in the analysis.
    """

    def __init__(
        self,
        intensities,
        lattice_ids,
        weights=None,
        min_pairs=3,
        lattice_group=None,
        dimensions=None,
        nproc=1,
        cc_weights=None,
    ):
        r"""Initialise a Target object.

        Args:
          intensities (cctbx.miller.array): The intensities on which to perform
            cosym analysis.
          lattice_ids (np.ndarray): An array of equal size to
            `intensities` which maps each reflection to a given lattice (dataset).
          weights (str): Optionally include weights in the target function.
            Allowed values are `None`, "count" and "standard_error". The default
            is to use no weights. If "count" is set, then weights are equal to the
            number of pairs of reflections used in calculating each value of the
            rij matrix. If "standard_error" is used, then weights are defined as
            :math:`w_{ij} = 1/s`, where :math:`s = (1-r_{ij}^2)/sqrt(N)`.
            Where N=(n-2) or N=(neff-1) depending on the cc_weights option.
            See also  https://doi.org/10.1525/collabra.87615.
          min_pairs (int): Only calculate the correlation coefficient between two
            datasets if they have more than `min_pairs` of common reflections.
          lattice_group (cctbx.sgtbx.space_group): Optionally set the lattice
            group to be used in the analysis.
          dimensions (int): Optionally override the number of dimensions to be used
            in the analysis. If not set, then the number of dimensions used is
            equal to the greater of 2 or the number of symmetry operations in the
            lattice group.
        """
        if weights is not None:
            assert weights in ("count", "standard_error")
        self._weights = weights
        self._min_pairs = min_pairs
        self._nproc = nproc

        data = intensities.customized_copy(anomalous_flag=False)
        cb_op_to_primitive = data.change_of_basis_op_to_primitive_setting()
        data = data.change_basis(cb_op_to_primitive).map_to_asu()

        # Convert to uint64 avoids crashes on Windows when later constructing
        # flex.size_t (https://github.com/cctbx/cctbx_project/issues/591)
        order = lattice_ids.argsort().astype(np.uint64)
        sorted_data = data.data().select(flex.size_t(order))
        sorted_indices = data.indices().select(flex.size_t(order))
        sorted_sigmas = data.sigmas().select(flex.size_t(order))
        self._lattice_ids = lattice_ids[order]
        self._data = data.customized_copy(
            indices=sorted_indices, data=sorted_data, sigmas=sorted_sigmas
        )
        assert isinstance(self._data.indices(), type(flex.miller_index()))
        assert isinstance(self._data.data(), type(flex.double()))

        # construct a lookup for the separate lattices
        self._lattices = np.array(
            [
                np.where(self._lattice_ids == i)[0][0]
                for i in np.unique(self._lattice_ids)
            ]
        )

        self.sym_ops = OrderedSet(["x,y,z"])
        self._lattice_group = lattice_group
        self.sym_ops.update(op.as_xyz() for op in self._generate_twin_operators())
        if dimensions is None:
            dimensions = max(2, len(self.sym_ops))
        self.set_dimensions(dimensions)

        self._lattice_group = copy.deepcopy(self._data.space_group())
        for sym_op in self.sym_ops:
            self._lattice_group.expand_smx(sym_op)
        self._patterson_group = self._lattice_group.build_derived_patterson_group()

        logger.debug(
            "Lattice group: %s (%i symops)",
            self._lattice_group.info().symbol_and_number(),
            len(self._lattice_group),
        )
        logger.debug(
            "Patterson group: %s", self._patterson_group.info().symbol_and_number()
        )
        if cc_weights == "sigma":
            self.rij_matrix, self.wij_matrix = self._compute_rij_wij_blockwise(
                cc_weights=True
            )
        else:
            self.rij_matrix, self.wij_matrix = self._compute_rij_wij_blockwise(
                cc_weights=False
            )

    def set_dimensions(self, dimensions):
        """Set the number of dimensions for analysis.

        Args:
          dimensions (int): The number of dimensions to be used.
        """
        self.dim = dimensions

    def _generate_twin_operators(self, lattice_symmetry_max_delta=5.0):
        # see also mmtbx.scaling.twin_analyses.twin_laws
        if self._lattice_group is None:
            cb_op_to_minimum_cell = self._data.change_of_basis_op_to_minimum_cell()
            minimum_cell_symmetry = self._data.crystal_symmetry().change_basis(
                cb_op=cb_op_to_minimum_cell
            )
            self._lattice_group = sgtbx.lattice_symmetry.group(
                reduced_cell=minimum_cell_symmetry.unit_cell(),
                max_delta=lattice_symmetry_max_delta,
            )
            intensity_symmetry = minimum_cell_symmetry.reflection_intensity_symmetry(
                anomalous_flag=self._data.anomalous_flag()
            )
            cb_op = cb_op_to_minimum_cell.inverse()
        else:
            cb_op = sgtbx.change_of_basis_op()
            intensity_symmetry = self._data.reflection_intensity_symmetry()

        operators = []
        for partition in cctbx.sgtbx.cosets.left_decomposition(
            g=self._lattice_group,
            h=intensity_symmetry.space_group()
            .build_derived_acentric_group()
            .make_tidy(),
        ).partitions[1:]:
            if partition[0].r().determinant() > 0:
                operators.append(cb_op.apply(partition[0]))

        return operators

    def _compute_rij_wij_blockwise(self, cc_weights=True):
        rij_matrix = None
        wij_matrix = None
        n = 0
        logger.info(
            f"Calculating rij matrix elements in {len(self._lattices)} row-blocks"
        )
        with concurrent.futures.ProcessPoolExecutor(max_workers=self._nproc) as pool:
            # note we use weights=True to help us work out where we have calculated rij,
            # even if the weights phil option is None
            futures = [
                pool.submit(
                    _compute_rij_matrix_one_row_block,
                    i,
                    self._lattices,
                    self._data,
                    self.sym_ops,
                    self._patterson_group,
                    weights=cc_weights,
                    min_pairs=self._min_pairs,
                )
                for i, _ in enumerate(self._lattices)
            ]
            for future in concurrent.futures.as_completed(futures):
                rij, wij = future.result()
                n += 1
                logger.info(f"Calculated rij matrix for row-block {n}")
                if rij_matrix is None:
                    rij_matrix = rij
                else:
                    rij_matrix += rij
                if wij is not None:
                    if wij_matrix is None:
                        wij_matrix = wij
                    else:
                        wij_matrix += wij

        rij_matrix = rij_matrix.toarray().astype(np.float64)
        rij_matrix += rij_matrix.T
        wij_matrix = wij_matrix.toarray().astype(np.float64)
        wij_matrix += wij_matrix.T

        ## Check if we have a dataset where no correlations could be calculate.
        zero_rows = np.where(np.all(wij_matrix == 0, axis=1))[0]
        if zero_rows.any():
            # only a problem if zero for all rows for dataset i
            Nlatt = len(self._lattices)
            Nsym = len(self.sym_ops)
            for i in range(Nlatt):
                all_rows = [i + j * Nlatt for j in range(Nsym)]
                if all(k in zero_rows for k in all_rows):
                    i_lower, i_upper = _lattice_lower_upper_index(self._lattices, i)
                    if not i_upper:
                        i_upper = len(self._data.data())
                    n = i_upper - i_lower
                    logger.warning(
                        f"Unable to calculate any correlations for datasets with index {i} ({n} reflections)."
                        + "\nIncreasing min_reflections or the resolution limit may overcome this problem."
                    )

        if self._weights:
            if self._weights == "standard_error":
                # N.B. using effective n due to sigma weighting, which can be below 2
                # but approches 1 in the limit, so rather say efective sample size
                # for standard error calc is n-1
                sel = np.where(wij_matrix > 1)
                se = (1 - np.square(rij_matrix[sel])) / np.sqrt(wij_matrix[sel] - 1)
                wij_matrix = np.zeros_like(rij_matrix)
                wij_matrix[sel] = 1 / se
            ## else uses the counts as weights
            # rescale the weights matrix such that the sum of wij_matrix == the number of non-zero entries
            scale = np.count_nonzero(wij_matrix) / np.sum(wij_matrix)
            wij_matrix *= scale
        else:
            ## No weights - i.e. equal weights in places where we can calculate an rij value,
            ## but also making sure our diagonal elements are zero as we exclude the
            ## self-correlation elements from rij and the cosym procedure - we need zero weights
            ## for uncalculate correlations so they aren't taken into account in the functional
            ## evaluation.
            ## at this point, wij matrix contains neff values where it was possible to calculate
            ## a pairwise correlation.
            sel = np.where(wij_matrix > 0)
            wij_matrix[sel] = 1

        return rij_matrix, wij_matrix

    def compute_functional(self, x: np.ndarray) -> float:
        """Compute the target function at coordinates `x`.

        Args:
          x (np.ndarray):
            a flattened list of the N-dimensional vectors, i.e. coordinates in
            the first dimension are stored first, followed by the coordinates in
            the second dimension, etc.

        Returns:
          f (float): The value of the target function at coordinates `x`.
        """
        assert (x.size // self.dim) == (len(self._lattices) * len(self.sym_ops))
        x = x.reshape((self.dim, x.size // self.dim))
        elements = np.square(self.rij_matrix - x.T @ x)
        if self.wij_matrix is not None:
            np.multiply(self.wij_matrix, elements, out=elements)
        f = 0.5 * elements.sum()
        return f

    def compute_functional_score_for_dimension_assessment(
        self, x: np.ndarray, outlier_rejection: bool = True
    ) -> float:
        if not outlier_rejection:
            return self.compute_functional(x)
        x = x.reshape((self.dim, x.size // self.dim))
        elements = np.square(self.rij_matrix - x.T @ x)
        if self.wij_matrix is not None:
            np.multiply(self.wij_matrix, elements, out=elements)

        q1, q2, q3 = np.quantile(elements, (0.25, 0.5, 0.75))
        inliers = elements[elements < q2 + (q3 - q1)]
        return 0.5 * inliers.sum()

    def compute_gradients_fd(self, x: np.ndarray, eps=1e-6) -> np.ndarray:
        """Compute the gradients at coordinates `x` using finite differences.

        Args:
          x (np.ndarray):
            a flattened list of the N-dimensional vectors, i.e. coordinates in
            the first dimension are stored first, followed by the coordinates in
            the second dimension, etc.
          eps (float):
            The value of epsilon to use in finite difference calculations.

        Returns:
          grad (np.ndarray):
          The gradients of the target function with respect to the parameters.
        """
        x = copy.deepcopy(x)
        grad = np.zeros(x.shape)
        for i in range(x.size):
            x[i] += eps  # x + eps
            fp = self.compute_functional(x)
            x[i] -= 2 * eps  # x - eps
            fm = self.compute_functional(x)
            x[i] += eps  # reset to original values
            grad[i] += (fp - fm) / (2 * eps)
        return grad

    def compute_gradients(self, x: np.ndarray) -> np.ndarray:
        """Compute the gradients of the target function at coordinates `x`.

        Args:
          x (np.ndarray):
            a flattened list of the N-dimensional vectors, i.e. coordinates in
            the first dimension are stored first, followed by the coordinates in
            the second dimension, etc.

        Returns:
          Tuple[float, np.ndarray]:
          f: The value of the target function at coordinates `x`.
          grad: The gradients of the target function with respect to the parameters.
        """
        x = x.reshape((self.dim, x.size // self.dim))
        if self.wij_matrix is not None:
            wrij_matrix = np.multiply(self.wij_matrix, self.rij_matrix)
            grad = -2 * x @ (wrij_matrix - np.multiply(self.wij_matrix, x.T @ x))
        else:
            grad = -2 * x @ (self.rij_matrix - x.T @ x)
        return grad.flatten()

    def curvatures(self, x: np.ndarray) -> np.ndarray:
        """Compute the curvature of the target function at coordinates `x`.

        Args:
          x (np.ndarray):
            a flattened list of the N-dimensional vectors, i.e. coordinates in
            the first dimension are stored first, followed by the coordinates in
            the second dimension, etc.

        Returns:
          curvs (np.ndarray):
          The curvature of the target function with respect to the parameters.
        """
        if self.wij_matrix is not None:
            wij = self.wij_matrix
        else:
            wij = np.ones(self.rij_matrix.shape)
        x = x.reshape((self.dim, x.size // self.dim))
        curvs = 2 * np.square(x) @ wij
        return curvs.flatten()

    def curvatures_fd(self, x: np.ndarray, eps=1e-6) -> np.ndarray:
        """Compute the curvatures at coordinates `x` using finite differences.

        Args:
          x (np.ndarray):
            a flattened list of the N-dimensional vectors, i.e. coordinates in
            the first dimension are stored first, followed by the coordinates in
            the second dimension, etc.
          eps (float):
            The value of epsilon to use in finite difference calculations.

        Returns:
          curvs (np.ndarray):
          The curvature of the target function with respect to the parameters.
        """
        x = copy.deepcopy(x)
        f = self.compute_functional(x)
        curvs = np.zeros(x.shape)
        for i in range(x.size):
            x[i] += eps  # x + eps
            fp = self.compute_functional(x)
            x[i] -= 2 * eps  # x - eps
            fm = self.compute_functional(x)
            x[i] += eps  # reset to original values
            curvs[i] += (fm - 2 * f + fp) / (eps**2)
        return curvs
