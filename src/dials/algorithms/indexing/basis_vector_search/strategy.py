"""Basis vector search strategies."""


from __future__ import annotations


class Strategy:
    """A base class for basis vector search strategies."""

    phil_help = None

    phil_scope = None

    def __init__(self, max_cell, params=None, *args, **kwargs):
        """Construct the strategy.

        Args:
            max_cell (float): An estimate of the maximum cell dimension of the primitive
                cell.
        """
        self._max_cell = max_cell
        self._params = params
        if self._params is None and self.phil_scope is not None:
            self._params = self.phil_scope.extract()

    def find_basis_vectors(self, reciprocal_lattice_vectors):
        """Find a list of likely basis vectors.

        Args:
            reciprocal_lattice_vectors (scitbx.array_family.flex.vec3_double):
                The list of reciprocal lattice vectors to search for periodicity.

        Returns:
            A tuple containing the list of basis vectors and a flex.bool array
            identifying which reflections were used in indexing.
        """
        raise NotImplementedError()
