from __future__ import annotations

import logging

import numpy

import iotbx.phil
from cctbx.sgtbx import space_group
from dxtbx import flumpy
from dxtbx.model import Crystal

from dials.algorithms.indexing import DialsIndexError

from .strategy import Strategy

# Import fast feedback indexer package (CUDA implementation of the TORO algorithm)
# https://github.com/paulscherrerinstitute/fast-feedback-indexer/tree/main/python
try:
    import ffbidx
except ModuleNotFoundError:
    ffbidx = None


logger = logging.getLogger(__name__)

toro_indexer_phil_str = """
toro_indexer
    .expert_level = 1
{
    max_output_cells = 32
        .type = int(value_min=1)
        .help = "Maximum number of output cells"
    max_spots = 300
        .type = int(value_min=8)
        .help = "Maximum number of reciprocal spots taken into account"
    num_candidate_vectors = 32
        .type = int(value_min=1)
        .help = "Number of candidate cell vectors"
    redundant_computations = True
        .type = bool
        .help = "Calculate candidates for all three cell vectors"
    dist1 = 0.3
        .type = float(value_min=0.001, value_max=0.5)
        .help = "Reciprocal spots within this threshold contribute to the score for vector sampling"
    dist3 = 0.15
        .type = float(value_min=0.001, value_max=0.8)
        .help = "Reciprocal spots within this threshold contribute to the score for cell sampling"
    num_halfsphere_points = 32768
        .type = int(value_min=8000)
        .help = "Number of sampling points on the half sphere"
    max_dist = 0.00075
        .type = float(value_min=0.0)
        .help = "Maximum final distance between measured and calculated reciprocal spot"
    min_spots = 8
        .type = int(value_min=6)
        .help = "Minimum number of reciprocal spots within distance max_dist"
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

        if ffbidx is None:
            raise DialsIndexError(
                "ToroIndexer requires fast feedback indexer. See (https://github.com/paulscherrerinstitute/fast-feedback-indexer)"
            )

        self._target_symmetry_primitive = target_symmetry_primitive
        self._max_lattices = max_lattices

        if target_symmetry_primitive is None:
            raise DialsIndexError("Target unit cell must be provided for TORO")

        target_cell = target_symmetry_primitive.unit_cell()
        if target_cell is None:
            raise ValueError("Please specify known_symmetry.unit_cell")

        self.params = params

        # Need the real space cell as numpy float32 array with all x vector coordinates, followed by y and z coordinates consecutively in memory
        self.input_cell = numpy.reshape(
            numpy.array(target_cell.fractionalization_matrix(), dtype="float32"), (3, 3)
        )

        # Create fast feedback indexer object (on default CUDA device)
        self.indexer = ffbidx.Indexer(
            max_output_cells=params.toro_indexer.max_output_cells,
            max_spots=params.toro_indexer.max_spots,
            num_candidate_vectors=params.toro_indexer.num_candidate_vectors,
            redundant_computations=params.toro_indexer.redundant_computations,
        )

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

        # Need the reciprocal lattice points as numpy float32 array with all x coordinates, followed by y and z coordinates consecutively in memory
        rlp = (
            numpy.array(flumpy.to_numpy(reflections["rlp"]), dtype="float32")
            .transpose()
            .copy()
        )

        output_cells, scores = self.indexer.run(
            rlp,
            self.input_cell,
            dist1=self.params.toro_indexer.dist1,
            dist3=self.params.toro_indexer.dist3,
            num_halfsphere_points=self.params.toro_indexer.num_halfsphere_points,
            max_dist=self.params.toro_indexer.max_dist,
            min_spots=self.params.toro_indexer.min_spots,
            n_output_cells=self.params.toro_indexer.max_output_cells,
        )

        cell_indices = self.indexer.crystals(
            output_cells,
            rlp,
            scores,
            threshold=self.params.toro_indexer.max_dist,
            min_spots=self.params.toro_indexer.min_spots,
        )

        candidate_crystal_models = []
        for index in cell_indices:
            j = 3 * index
            real_a = output_cells[:, j]
            real_b = output_cells[:, j + 1]
            real_c = output_cells[:, j + 2]
            crystal = Crystal(real_a, real_b, real_c, space_group=space_group("P1"))
            candidate_crystal_models.append(crystal)

        return candidate_crystal_models
