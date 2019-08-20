from __future__ import absolute_import, division, print_function

from libtbx import phil
from scitbx.array_family import flex
from scitbx import matrix

from . import strategies


fft1d_phil_str = """\
characteristic_grid = None
    .help = Sampling frequency in radians. See Steller 1997. If None, \
            determine a grid sampling automatically using the input \
            reflections, using at most 0.029 radians.
    .type = float(value_min=0)
"""


class FFT1D(strategies.Strategy):
    """Basis vector search using a 1D FFT.

    See:
        Steller, I., Bolotovsky, R. & Rossmann, M. G. (1997). J. Appl. Cryst. 30, 1036-1040.
        Sauter, N. K., Grosse-Kunstleve, R. W. & Adams, P. D. (2004). J. Appl. Cryst. 37, 399-409.

    """

    phil_scope = phil.parse(fft1d_phil_str)

    def __init__(self, max_cell, params=None, *args, **kwargs):
        """Construct an FFT1D object.

        Args:
            max_cell (float): An estimate of the maximum cell dimension of the primitive
                cell.
            characteristic_grid (float): Sampling frequency in radians. See Steller 1997.
                If None, determine a grid sampling automatically using the input
                reflections, using at most 0.029 radians.

        """
        super(FFT1D, self).__init__(max_cell, params=params, *args, **kwargs)

    def find_basis_vectors(self, reciprocal_lattice_vectors):
        """Find a list of likely basis vectors.

        Args:
            reciprocal_lattice_vectors (scitbx.array_family.flex.vec3_double):
                The list of reciprocal lattice vectors to search for periodicity.

        Returns:
            A tuple containing the list of basis vectors and a flex.bool array
            identifying which reflections were used in indexing.

        """
        from rstbx.phil.phil_preferences import indexing_api_defs
        import iotbx.phil

        used_in_indexing = flex.bool(reciprocal_lattice_vectors.size(), True)

        hardcoded_phil = iotbx.phil.parse(input_string=indexing_api_defs).extract()

        # Spot_positions: Centroid positions for spotfinder spots, in pixels
        # Return value: Corrected for parallax, converted to mm

        # derive a max_cell from mm spots
        # derive a grid sampling from spots

        from rstbx.indexing_api.lattice import DPS_primitive_lattice

        # max_cell: max possible cell in Angstroms; set to None, determine from data
        # recommended_grid_sampling_rad: grid sampling in radians; guess for now

        DPS = DPS_primitive_lattice(
            max_cell=self._max_cell,
            recommended_grid_sampling_rad=self._params.characteristic_grid,
            horizon_phil=hardcoded_phil,
        )

        # transform input into what Nick needs
        # i.e., construct a flex.vec3 double consisting of mm spots, phi in degrees
        DPS.index(reciprocal_space_vectors=reciprocal_lattice_vectors)
        solutions = DPS.getSolutions()
        candidate_basis_vectors = [matrix.col(s.bvec()) for s in solutions]
        return candidate_basis_vectors, used_in_indexing
