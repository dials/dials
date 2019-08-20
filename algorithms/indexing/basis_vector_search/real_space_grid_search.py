from __future__ import absolute_import, division, print_function

import logging
import math

from libtbx import phil
from scitbx import matrix

from . import strategies


logger = logging.getLogger(__name__)


real_space_grid_search_phil_str = """\
characteristic_grid = 0.02
    .type = float(value_min=0)
"""


class RealSpaceGridSearch(strategies.Strategy):
    """Basis vector search using a real space grid search.

    See:
        Gildea, R. J., Waterman, D. G., Parkhurst, J. M., Axford, D., Sutton, G., Stuart, D. I., Sauter, N. K., Evans, G. & Winter, G. (2014). Acta Cryst. D70, 2652-2666.

    """

    phil_scope = phil.parse(real_space_grid_search_phil_str)

    def __init__(self, max_cell, target_unit_cell, params=None, *args, **kwargs):
        """Construct a real_space_grid_search object.

        Args:
            max_cell (float): An estimate of the maximum cell dimension of the primitive
                cell.
            target_unit_cell (cctbx.uctbx.unit_cell): The target unit cell.
            characteristic_grid (float): Sampling frequency in radians.

        """
        super(RealSpaceGridSearch, self).__init__(
            max_cell, params=params, *args, **kwargs
        )
        self._target_unit_cell = target_unit_cell

    def find_basis_vectors(self, reciprocal_lattice_vectors):
        """Find a list of likely basis vectors.

        Args:
            reciprocal_lattice_vectors (scitbx.array_family.flex.vec3_double):
                The list of reciprocal lattice vectors to search for periodicity.

        """
        from rstbx.array_family import (
            flex,
        )  # required to load scitbx::af::shared<rstbx::Direction> to_python converter

        used_in_indexing = flex.bool(reciprocal_lattice_vectors.size(), True)

        logger.info("Indexing from %i reflections" % used_in_indexing.count(True))

        def compute_functional(vector):
            two_pi_S_dot_v = 2 * math.pi * reciprocal_lattice_vectors.dot(vector)
            return flex.sum(flex.cos(two_pi_S_dot_v))

        from rstbx.dps_core import SimpleSamplerTool

        SST = SimpleSamplerTool(self._params.characteristic_grid)
        SST.construct_hemisphere_grid(SST.incr)
        cell_dimensions = self._target_unit_cell.parameters()[:3]
        unique_cell_dimensions = set(cell_dimensions)
        logger.info(
            "Number of search vectors: %i"
            % (len(SST.angles) * len(unique_cell_dimensions))
        )
        vectors = flex.vec3_double()
        function_values = flex.double()
        for i, direction in enumerate(SST.angles):
            for l in unique_cell_dimensions:
                v = matrix.col(direction.dvec) * l
                f = compute_functional(v.elems)
                vectors.append(v.elems)
                function_values.append(f)

        perm = flex.sort_permutation(function_values, reverse=True)
        vectors = vectors.select(perm)
        function_values = function_values.select(perm)

        unique_vectors = []
        i = 0
        while len(unique_vectors) < 30:
            v = matrix.col(vectors[i])
            is_unique = True
            if i > 0:
                for v_u in unique_vectors:
                    if v.length() < v_u.length():
                        if strategies._is_approximate_integer_multiple(v, v_u):
                            is_unique = False
                            break
                    elif strategies._is_approximate_integer_multiple(v_u, v):
                        is_unique = False
                        break
            if is_unique:
                unique_vectors.append(v)
            i += 1

        for i in range(30):
            v = matrix.col(vectors[i])
            logger.debug(
                "%s %s %s" % (str(v.elems), str(v.length()), str(function_values[i]))
            )

        logger.info("Number of unique vectors: %i" % len(unique_vectors))

        for i in range(len(unique_vectors)):
            logger.debug(
                "%s %s %s"
                % (
                    str(compute_functional(unique_vectors[i].elems)),
                    str(unique_vectors[i].length()),
                    str(unique_vectors[i].elems),
                )
            )

        return unique_vectors, used_in_indexing
