from __future__ import annotations

import logging
import math

from libtbx import phil
from rstbx.array_family import (
    flex,  # required to load scitbx::af::shared<rstbx::Direction> to_python converter
)
from rstbx.dps_core import SimpleSamplerTool
from scitbx import matrix

from dials.algorithms.indexing import DialsIndexError

from .strategy import Strategy
from .utils import group_vectors

logger = logging.getLogger(__name__)


real_space_grid_search_phil_str = """\
characteristic_grid = 0.02
    .type = float(value_min=0)
max_vectors = 30
    .help = "The maximum number of unique vectors to find in the grid search."
    .type = int(value_min=3)
"""


class RealSpaceGridSearch(Strategy):
    """
    Basis vector search using a real space grid search.

    Search strategy to index found spots based on known unit cell parameters. It is
    often useful for difficult cases of narrow-wedge rotation data or stills data,
    especially where there is diffraction from multiple crystals.

    A set of dimensionless radial unit vectors, typically ~7000 in total, is chosen
    so that they are roughly evenly spaced in solid angle over a hemisphere. For each
    direction, each of the three known unit cell vectors is aligned with the unit
    vector and is scored according to how well it accords with the periodicity in
    that direction of the reconstructed reciprocal space positions of the observed
    spot centroids. Examining the highest-scoring combinations, any basis vectors in
    orientations that are nearly collinear with a shorter basis vector are
    eliminated. The highest-scoring remaining combinations are selected as the basis
    of the direct lattice.

    See:
        Gildea, R. J., Waterman, D. G., Parkhurst, J. M., Axford, D., Sutton, G., Stuart, D. I., Sauter, N. K., Evans, G. & Winter, G. (2014). Acta Cryst. D70, 2652-2666.
    """

    phil_help = (
        "Index the found spots by testing a known unit cell in various orientations "
        "until the best match is found. This strategy is often useful for difficult "
        "cases of narrow-wedge rotation data or stills data, especially where there "
        "is diffraction from multiple crystals."
    )

    phil_scope = phil.parse(real_space_grid_search_phil_str)

    def __init__(self, max_cell, target_unit_cell, params=None, *args, **kwargs):
        """Construct a real_space_grid_search object.

        Args:
            max_cell (float): An estimate of the maximum cell dimension of the primitive
                cell.
            target_unit_cell (cctbx.uctbx.unit_cell): The target unit cell.
        """
        super().__init__(max_cell, params=params, *args, **kwargs)
        if target_unit_cell is None:
            raise DialsIndexError(
                "Target unit cell must be provided for real_space_grid_search"
            )
        self._target_unit_cell = target_unit_cell

    @property
    def search_directions(self):
        """Generator of the search directions (i.e. vectors with length 1)."""
        SST = SimpleSamplerTool(self._params.characteristic_grid)
        SST.construct_hemisphere_grid(SST.incr)
        for direction in SST.angles:
            yield matrix.col(direction.dvec)

    @property
    def search_vectors(self):
        """Generator of the search vectors.

        The lengths of the vectors correspond to the target unit cell dimensions.
        """
        unique_cell_dimensions = set(self._target_unit_cell.parameters()[:3])
        for i, direction in enumerate(self.search_directions):
            for l in unique_cell_dimensions:
                yield direction * l

    @staticmethod
    def compute_functional(vector, reciprocal_lattice_vectors):
        """Compute the functional for a single direction vector.

        Args:
            vector (tuple): The vector at which to compute the functional.
            reciprocal_lattice_vectors (scitbx.array_family.flex.vec3_double):
                The list of reciprocal lattice vectors.

        Returns:
            The functional for the given vector.
        """
        two_pi_S_dot_v = 2 * math.pi * reciprocal_lattice_vectors.dot(vector)
        return flex.sum(flex.cos(two_pi_S_dot_v))

    def score_vectors(self, reciprocal_lattice_vectors):
        """Compute the functional for the given directions.

        Args:
            directions: An iterable of the search directions.
            reciprocal_lattice_vectors (scitbx.array_family.flex.vec3_double):
                The list of reciprocal lattice vectors.
        Returns:
            A tuple containing the list of search vectors and their scores.
        """
        vectors = flex.vec3_double()
        scores = flex.double()
        for i, v in enumerate(self.search_vectors):
            f = self.compute_functional(v.elems, reciprocal_lattice_vectors)
            vectors.append(v.elems)
            scores.append(f)
        return vectors, scores

    def find_basis_vectors(self, reciprocal_lattice_vectors):
        """Find a list of likely basis vectors.

        Args:
            reciprocal_lattice_vectors (scitbx.array_family.flex.vec3_double):
                The list of reciprocal lattice vectors to search for periodicity.
        """
        used_in_indexing = flex.bool(reciprocal_lattice_vectors.size(), True)
        logger.info("Indexing from %i reflections", used_in_indexing.count(True))

        vectors, weights = self.score_vectors(reciprocal_lattice_vectors)

        perm = flex.sort_permutation(weights, reverse=True)
        vectors = vectors.select(perm)
        weights = weights.select(perm)

        groups = group_vectors(vectors, weights, max_groups=self._params.max_vectors)
        unique_vectors = []
        unique_weights = []
        for g in groups:
            idx = flex.max_index(flex.double(g.weights))
            unique_vectors.append(g.vectors[idx])
            unique_weights.append(g.weights[idx])

        logger.info("Number of unique vectors: %i", len(unique_vectors))

        for v, w in zip(unique_vectors, unique_weights):
            logger.debug("%s %s %s", w, v.length(), str(v.elems))

        return unique_vectors, used_in_indexing
