from __future__ import annotations

import logging
import math

from cctbx import crystal, uctbx, xray
from libtbx import libtbx, phil
from scitbx import fftpack, matrix
from scitbx.array_family import flex

import dials_algorithms_indexing_ext
from dials.algorithms import indexing

from .strategy import Strategy
from .utils import group_vectors, is_approximate_integer_multiple

logger = logging.getLogger(__name__)


fft3d_phil_str = """\
b_iso = Auto
    .type = float(value_min=0)
    .expert_level = 2
rmsd_cutoff = 15
    .type = float(value_min=0)
    .expert_level = 1
peak_search = *flood_fill clean
    .type = choice
    .expert_level = 2
peak_volume_cutoff = 0.15
    .type = float
    .expert_level = 2
reciprocal_space_grid {
    n_points = 256
        .type = int(value_min=0)
        .expert_level = 1
    d_min = Auto
        .type = float(value_min=0)
        .help = "The high resolution limit in Angstrom for spots to include in "
                "the initial indexing."
    }
"""


class FFT3D(Strategy):
    """
    Basis vector search using a 3D FFT in reciprocal space.

    Reciprocal space is sampled as a 3D Cartesian grid, aligned with the basis of the
    laboratory frame. The reciprocal-space positions of the centroids of measured
    spots are ascribed a value of 1, with the rest of the grid assigned a value of 0.
    A 3D FFT is performed and the three shortest non-collinear reciprocal spatial
    wave vectors with appreciable spectral weight correspond to the basis vectors of
    the real space lattice.

    Because this procedure requires a sampling of all of reciprocal space, up to the
    d* value of the measured spot with the highest resolution, it can be more memory
    intensive than alternative approaches. To mitigate this, the 3D FFT will
    sometimes be curtailed to a region of reciprocal space below a certain
    resolution, and higher-resolution spots will be ignored.

    See:
        Bricogne, G. (1986). Proceedings of the EEC Cooperative Workshop on Position-Sensitive Detector Software (Phase III), p. 28. Paris: LURE.
        Campbell, J. W. (1998). J. Appl. Cryst. 31, 407-413.
    """

    phil_help = (
        "Search for the basis vectors of the direct lattice by performing a 3D FFT in "
        "reciprocal space of the density of found spots. Since this can be quite "
        "memory-intensive, the data used for indexing may automatically be "
        "constrained to just the lower resolution spots."
    )

    phil_scope = phil.parse(fft3d_phil_str)

    def __init__(self, max_cell, min_cell=3, params=None, *args, **kwargs):
        """Construct an FFT3D object.

        Args:
            max_cell (float): An estimate of the maximum cell dimension of the primitive
                cell.
            n_points (int): The size of the fft3d grid.
            d_min (float): The high resolution limit in Angstrom for spots to include in
                the initial indexing. If `Auto` then calculated as
                `d_min = 5 * max_cell/n_points`.
            b_iso (float): Apply an isotropic b_factor weight to the points when doing
                the FFT. If `Auto` then calculated as
                `b_iso = -4 * d_min ** 2 * math.log(0.05)`.
            rmsd_cutoff (float): RMSD cutoff applied to the transformed real-space map
                prior to performing the peak search.
            peak_volume_cutoff (float): Only include peaks that are larger than this
                fraction of the volume of the largest peak in the transformed real-space
                map.
            min_cell (float): A conservative lower bound on the minimum possible
                primitive unit cell dimension.
        """
        super().__init__(max_cell, params=params, *args, **kwargs)
        n_points = self._params.reciprocal_space_grid.n_points
        self._gridding = fftpack.adjust_gridding_triple(
            (n_points, n_points, n_points), max_prime=5
        )
        self._n_points = self._gridding[0]
        self._min_cell = min_cell

    def find_basis_vectors(self, reciprocal_lattice_vectors):
        """Find a list of likely basis vectors.

        Args:
            reciprocal_lattice_vectors (scitbx.array_family.flex.vec3_double):
                The list of reciprocal lattice vectors to search for periodicity.

        Returns:
            A tuple containing the list of basis vectors and a flex.bool array
            identifying which reflections were used in indexing.
        """
        if self._params.reciprocal_space_grid.d_min is libtbx.Auto:
            # rough calculation of suitable d_min based on max cell
            # see also Campbell, J. (1998). J. Appl. Cryst., 31(3), 407-413.
            # fft_cell should be greater than twice max_cell, so say:
            #   fft_cell = 2.5 * max_cell
            # then:
            #   fft_cell = n_points * d_min/2
            #   2.5 * max_cell = n_points * d_min/2
            # a little bit of rearrangement:
            #   d_min = 5 * max_cell/n_points

            max_cell = self._max_cell
            d_min = 5 * max_cell / self._n_points
            d_spacings = 1 / reciprocal_lattice_vectors.norms()
            d_min = max(d_min, min(d_spacings))
            logger.info("Setting d_min: %.2f", d_min)
        else:
            d_min = self._params.reciprocal_space_grid.d_min

        grid_real, used_in_indexing = self._fft(reciprocal_lattice_vectors, d_min)
        self.sites, self.volumes = self._find_peaks(grid_real, d_min)

        # hijack the xray.structure class to facilitate calculation of distances

        self.crystal_symmetry = crystal.symmetry(
            unit_cell=self._fft_cell, space_group_symbol="P1"
        )
        xs = xray.structure(crystal_symmetry=self.crystal_symmetry)
        for i, site in enumerate(self.sites):
            xs.add_scatterer(xray.scatterer("C%i" % i, site=site))

        xs = xs.sites_mod_short()
        sites_cart = xs.sites_cart()
        lengths = flex.double([matrix.col(sc).length() for sc in sites_cart])
        perm = flex.sort_permutation(lengths)
        xs = xs.select(perm)
        volumes = self.volumes.select(perm)

        vectors = xs.sites_cart()
        norms = vectors.norms()
        sel = (norms > self._min_cell) & (norms < (2 * self._max_cell))
        vectors = vectors.select(sel)
        vectors = [matrix.col(v) for v in vectors]
        volumes = volumes.select(sel)

        vector_groups = group_vectors(vectors, volumes)
        vectors = [g.mean for g in vector_groups]
        volumes = flex.double(max(g.weights) for g in vector_groups)

        # sort by peak size
        perm = flex.sort_permutation(volumes, reverse=True)
        volumes = volumes.select(perm)
        vectors = [vectors[i] for i in perm]

        for i, (v, volume) in enumerate(zip(vectors, volumes)):
            logger.debug(f"{i} {v.length()} {volume}")

        # sort by length
        lengths = flex.double(v.length() for v in vectors)
        perm = flex.sort_permutation(lengths)

        # exclude vectors that are (approximately) integer multiples of a shorter
        # vector
        unique_vectors = []
        unique_volumes = flex.double()
        for p in perm:
            v = vectors[p]
            is_unique = True
            for i, v_u in enumerate(unique_vectors):
                if (unique_volumes[i] > volumes[p]) and is_approximate_integer_multiple(
                    v_u, v
                ):
                    logger.debug(
                        "rejecting %s: integer multiple of %s", v.length(), v_u.length()
                    )
                    is_unique = False
                    break
            if is_unique:
                unique_vectors.append(v)
                unique_volumes.append(volumes[p])

        # re-sort by peak volume
        perm = flex.sort_permutation(unique_volumes, reverse=True)
        self.candidate_basis_vectors = [unique_vectors[i] for i in perm]
        return self.candidate_basis_vectors, used_in_indexing

    def _fft(self, reciprocal_lattice_vectors, d_min):

        (
            reciprocal_space_grid,
            used_in_indexing,
        ) = self._map_centroids_to_reciprocal_space_grid(
            reciprocal_lattice_vectors, d_min
        )

        logger.info(
            "Number of centroids used: %i", (reciprocal_space_grid > 0).count(True)
        )

        # gb_to_bytes = 1073741824
        # bytes_to_gb = 1/gb_to_bytes
        # (128**3)*8*2*bytes_to_gb
        # 0.03125
        # (256**3)*8*2*bytes_to_gb
        # 0.25
        # (512**3)*8*2*bytes_to_gb
        # 2.0

        fft = fftpack.complex_to_complex_3d(self._gridding)
        grid_complex = flex.complex_double(
            reals=reciprocal_space_grid,
            imags=flex.double(reciprocal_space_grid.size(), 0),
        )
        grid_transformed = fft.forward(grid_complex)
        grid_real = flex.pow2(flex.real(grid_transformed))
        del grid_transformed

        return grid_real, used_in_indexing

    def _map_centroids_to_reciprocal_space_grid(
        self, reciprocal_lattice_vectors, d_min
    ):
        logger.info("FFT gridding: (%i,%i,%i)" % self._gridding)

        grid = flex.double(flex.grid(self._gridding), 0)

        if self._params.b_iso is libtbx.Auto:
            self._params.b_iso = -4 * d_min**2 * math.log(0.05)
            logger.debug("Setting b_iso = %.1f", self._params.b_iso)
        used_in_indexing = flex.bool(reciprocal_lattice_vectors.size(), True)
        dials_algorithms_indexing_ext.map_centroids_to_reciprocal_space_grid(
            grid,
            reciprocal_lattice_vectors,
            used_in_indexing,  # do we really need this?
            d_min,
            b_iso=self._params.b_iso,
        )
        return grid, used_in_indexing

    def _find_peaks(self, grid_real, d_min):
        grid_real_binary = grid_real.deep_copy()
        rmsd = math.sqrt(
            flex.mean(
                flex.pow2(
                    grid_real_binary.as_1d() - flex.mean(grid_real_binary.as_1d())
                )
            )
        )
        grid_real_binary.set_selected(
            grid_real_binary < (self._params.rmsd_cutoff) * rmsd, 0
        )
        grid_real_binary.as_1d().set_selected(grid_real_binary.as_1d() > 0, 1)
        grid_real_binary = grid_real_binary.iround()
        from cctbx import masks

        # real space FFT grid dimensions
        cell_lengths = [self._n_points * d_min / 2 for i in range(3)]
        self._fft_cell = uctbx.unit_cell(cell_lengths + [90] * 3)

        flood_fill = masks.flood_fill(grid_real_binary, self._fft_cell)
        if flood_fill.n_voids() < 4:
            # Require at least peak at origin and one peak for each basis vector
            raise indexing.DialsIndexError(
                "Indexing failed: fft3d peak search failed to find sufficient number of peaks."
            )

        # the peak at the origin might have a significantly larger volume than the
        # rest so exclude any anomalously large peaks from determining minimum volume
        from scitbx.math import five_number_summary

        outliers = flex.bool(flood_fill.n_voids(), False)
        grid_points_per_void = flood_fill.grid_points_per_void()
        min_x, q1_x, med_x, q3_x, max_x = five_number_summary(grid_points_per_void)
        iqr_multiplier = 5
        iqr_x = q3_x - q1_x
        cut_x = iqr_multiplier * iqr_x
        outliers.set_selected(grid_points_per_void.as_double() > (q3_x + cut_x), True)
        # print q3_x + cut_x, outliers.count(True)
        isel = (
            grid_points_per_void
            > int(
                self._params.peak_volume_cutoff
                * flex.max(grid_points_per_void.select(~outliers))
            )
        ).iselection()

        sites = flood_fill.centres_of_mass_frac().select(isel)
        volumes = flood_fill.grid_points_per_void().select(isel)
        return sites, volumes
