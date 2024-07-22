from __future__ import annotations

import logging

import iotbx.phil

from dials.algorithms.indexing import DialsIndexError

from .strategy import Strategy

logger = logging.getLogger(__name__)

toro_indexer_phil_str = """
toro_indexer
    .expert_level = 1
{
    some_parameter = 1
        .type = int
}
"""


class ToroIndexer(Strategy):
    """
    A lattice search strategy using the TORO algorithm.
    For more info, see:
    [Gasparotto P, et al. TORO Indexer: a PyTorch-based indexing algorithm for kilohertz serial crystallography. J. Appl. Cryst. 2024 57(4)](https://doi.org/10.1107/S1600576724003182)
    """

    phil_help = (
        "A lattice search strategy for very fast indexing using PyTorch acceleration"
    )

    phil_scope = iotbx.phil.parse(toro_indexer_phil_str)

    def __init__(
        self, target_symmetry_primitive, max_lattices, params=None, *args, **kwargs
    ):
        """Construct ToroIndexer object.

        Args:
            target_symmetry_primitive (cctbx.crystal.symmetry): The target
                crystal symmetry and unit cell
            max_lattices (int): The maximum number of lattice models to find
            params (phil,optional): Phil params

        Returns:
            None
        """
        super().__init__(params=None, *args, **kwargs)
        self._target_symmetry_primitive = target_symmetry_primitive
        self._max_lattices = max_lattices

        if target_symmetry_primitive is None:
            raise DialsIndexError(
                "Target unit cell and space group must be provided for TORO"
            )

        target_cell = target_symmetry_primitive.unit_cell()
        if target_cell is None:
            raise ValueError("Please specify known_symmetry.unit_cell")

    def find_crystal_models(self, reflections, experiments):
        """Find a list of candidate crystal models.

        Args:
            reflections (dials.array_family.flex.reflection_table):
                The found spots centroids and associated data

            experiments (dxtbx.model.experiment_list.ExperimentList):
                The experimental geometry models

        Returns:
            A list of candidate crystal models.
        """

        raise NotImplementedError
