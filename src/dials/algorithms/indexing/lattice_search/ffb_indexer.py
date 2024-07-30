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

ffbidx_phil_str = """
ffbidx
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
    method = *ifssr ifss ifse raw
        .type = choice
        .help = "Refinement method (consult algorithm description)"
    triml = 0.001
        .type = float(value_min=0, value_max=0.5)
        .help = "lower trimming value for intermediate score calculations"
    trimh = 0.3
        .type = float(value_min=0, value_max=0.5)
        .help = "higher trimming value for intermediate score calculations"
    delta = 0.1
        .type = float(value_min=0.000001)
        .help = "log2 curve position for intermediate score calculations, lower values will me more selective in choosing close spots"
    simple_data_filename = None
        .type = path
        .help = "Optional filename for the output of a simple data file for debugging"
        .expert_level = 2
}
"""


def write_simple_data_file(filename, rlp, cell):
    """Write a simple data file for debugging."""
    with open(filename, "w") as f:
        f.write(" ".join(map(str, cell.ravel())) + "\n")
        for r in rlp:
            f.write(" ".join(map(str, r.ravel())) + "\n")


class FfbIndexer(Strategy):
    """
    A lattice search strategy using a Cuda-accelerated implementation of the TORO algorithm.
    For more info, see:
    [Gasparotto P, et al. TORO Indexer: a PyTorch-based indexing algorithm for kilohertz serial crystallography. J. Appl. Cryst. 2024 57(4)](https://doi.org/10.1107/S1600576724003182)
    """

    phil_help = (
        "A lattice search strategy for very fast indexing using Cuda acceleration"
    )

    phil_scope = iotbx.phil.parse(ffbidx_phil_str)

    def __init__(
        self, target_symmetry_primitive, max_lattices, params=None, *args, **kwargs
    ):
        """Construct FfbIndexer object.

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
                "ffbidx requires the fast feedback indexer package. See (https://github.com/paulscherrerinstitute/fast-feedback-indexer)"
            )

        self._target_symmetry_primitive = target_symmetry_primitive
        self._max_lattices = max_lattices

        if target_symmetry_primitive is None:
            raise DialsIndexError("Target unit cell must be provided for ffbidx")

        target_cell = target_symmetry_primitive.unit_cell()
        if target_cell is None:
            raise ValueError("Please specify known_symmetry.unit_cell")

        self.params = params

        # Need the real space cell as numpy float32 array with all x vector coordinates, followed by y and z coordinates consecutively in memory
        self.input_cell = numpy.reshape(
            numpy.array(target_cell.orthogonalization_matrix(), dtype="float32"), (3, 3)
        )

        # Create fast feedback indexer object (on default CUDA device)
        self.indexer = ffbidx.Indexer(
            max_output_cells=params.ffbidx.max_output_cells,
            max_spots=params.ffbidx.max_spots,
            num_candidate_vectors=params.ffbidx.num_candidate_vectors,
            redundant_computations=params.ffbidx.redundant_computations,
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

        rlp = numpy.array(flumpy.to_numpy(reflections["rlp"]), dtype="float32")

        if self.params.ffbidx.simple_data_filename is not None:
            write_simple_data_file(
                self.params.ffbidx.simple_data_filename, rlp, self.input_cell
            )

        # Need the reciprocal lattice points as numpy float32 array with all x coordinates, followed by y and z coordinates consecutively in memory
        rlp = rlp.transpose().copy()

        output_cells, scores = self.indexer.run(
            rlp,
            self.input_cell,
            dist1=self.params.ffbidx.dist1,
            dist3=self.params.ffbidx.dist3,
            num_halfsphere_points=self.params.ffbidx.num_halfsphere_points,
            max_dist=self.params.ffbidx.max_dist,
            min_spots=self.params.ffbidx.min_spots,
            n_output_cells=self.params.ffbidx.max_output_cells,
            method=self.params.ffbidx.method,
            triml=self.params.ffbidx.triml,
            trimh=self.params.ffbidx.trimh,
            delta=self.params.ffbidx.delta,
        )

        cell_indices = self.indexer.crystals(
            output_cells,
            rlp,
            scores,
            threshold=self.params.ffbidx.max_dist,
            min_spots=self.params.ffbidx.min_spots,
            method=self.params.ffbidx.method,
        )

        candidate_crystal_models = []
        if cell_indices is None:
            return candidate_crystal_models

        for index in cell_indices:
            j = 3 * index
            real_a = output_cells[:, j]
            real_b = output_cells[:, j + 1]
            real_c = output_cells[:, j + 2]
            crystal = Crystal(
                real_a.tolist(),
                real_b.tolist(),
                real_c.tolist(),
                space_group=space_group("P1"),
            )
            candidate_crystal_models.append(crystal)

        return candidate_crystal_models
