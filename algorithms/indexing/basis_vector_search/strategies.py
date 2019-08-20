"""Basis vector search strategies."""

from __future__ import absolute_import, division, print_function

import abc


class Strategy(object):
    """A base class for basis vector search strategies."""

    __metaclass__ = abc.ABCMeta

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

    @abc.abstractmethod
    def find_basis_vectors(self, reciprocal_lattice_vectors):
        """Find a list of likely basis vectors.

        Args:
            reciprocal_lattice_vectors (scitbx.array_family.flex.vec3_double):
                The list of reciprocal lattice vectors to search for periodicity.

        Returns:
            A tuple containing the list of basis vectors and a flex.bool array
            identifying which reflections were used in indexing.

        """
        pass


def _is_approximate_integer_multiple(
    vec_a, vec_b, relative_tolerance=0.2, angular_tolerance=5.0
):
    length_a = vec_a.length()
    length_b = vec_b.length()
    # assert length_b >= length_a
    if length_a > length_b:
        vec_a, vec_b = vec_b, vec_a
        length_a, length_b = length_b, length_a
    angle = vec_a.angle(vec_b, deg=True)
    if angle < angular_tolerance or abs(180 - angle) < angular_tolerance:
        n = length_b / length_a
        if abs(round(n) - n) < relative_tolerance:
            return True
    return False
