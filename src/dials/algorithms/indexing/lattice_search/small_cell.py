from __future__ import annotations

# import copy
import logging

# import math
# import operator
#
import libtbx.phil

# from dials.array_family import flex
from xfel.small_cell.small_cell import small_cell_index_detail

# from cctbx import miller
# from dxtbx.model import Crystal
# from scitbx import matrix
# from scitbx.math import least_squares_plane, superpose
#
from dials.algorithms.indexing import DialsIndexError

from .strategy import Strategy

logger = logging.getLogger(__name__)


small_cell_phil_str = """\
include scope xfel.small_cell.command_line.small_cell_index.small_cell_phil_str
"""


class SmallCell(Strategy):
    """
    A lattice search strategy using the small_cell algorithm.

    """

    phil_help = (
        "A lattice search strategy that matches low resolution spots to candidate "
        "indices based on a known unit cell and space group."
    )

    phil_scope = libtbx.phil.parse(small_cell_phil_str)

    def __init__(
        self, target_symmetry_primitive, max_lattices, params=None, *args, **kwargs
    ):
        """Construct a LowResSpotMatch object.

        Args:
            target_symmetry_primitive (cctbx.crystal.symmetry): The target
                crystal symmetry and unit cell

            max_lattices (int): The maximum number of lattice models to find
        """
        super().__init__(params=params, *args, **kwargs)
        self._target_symmetry_primitive = target_symmetry_primitive
        self._max_lattices = max_lattices

        if target_symmetry_primitive is None:
            raise DialsIndexError(
                "Target unit cell and space group must be provided for small_cell"
            )

        params.small_cell.powdercell = (
            target_symmetry_primitive.unit_cell().parameters()
        )
        params.small_cell.spacegroup = str(target_symmetry_primitive.space_group_info())

    def find_crystal_models(self, reflections, experiments):
        """Find a list of candidate crystal models.

        Args:
            reflections (dials.array_family.flex.reflection_table):
                The found spots centroids and associated data

            experiments (dxtbx.model.experiment_list.ExperimentList):
                The experimental geometry models
        """

        max_clique_len, indexed_experiments, refls = small_cell_index_detail(
            experiments, reflections, self._params, write_output=False
        )

        # FIXME small_cell_index_detail creates an indexed reflections list, but
        # we don't use that in dials.index...

        candidate_crystal_models = []

        self.candidate_crystal_models = candidate_crystal_models
        return self.candidate_crystal_models
