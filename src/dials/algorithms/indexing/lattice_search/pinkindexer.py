from __future__ import annotations

import logging

import gemmi
from scipy.linalg import orthogonal_procrustes
from scipy.optimize import minimize

import iotbx.phil
from dxtbx.model import Crystal
from scitbx import matrix

from dials.algorithms.indexing import DialsIndexError

from .strategy import Strategy

logger = logging.getLogger(__name__)

import numpy as np

pink_indexer_phil_str = """
pink_indexer
    .expert_level = 1
{
    max_refls = 50
        .type = int(value_min=10)
        .help = "Maximum number of reflections to consider indexing"
    wav_min = None
        .type = float(value_min=0.)
        .help = "Minimum wavelength for polychromatic data"
    wav_max = None
        .type = float(value_min=0.)
        .help = "Maximum wavelength for polychromatic data"
    max_refls = 50
        .type = int(value_min=10, value_max=100)
        .help = "Maximum number of reflections to consider during indexing"
}
"""


def angle_between(vec1, vec2, deg=True):
    """
    This function computes the angle between vectors along the last dimension of the input arrays.
    This version is a numerically stable one based on arctan2 as described in this post:
     - https://scicomp.stackexchange.com/a/27769/39858

    Parameters
    ----------
    vec1 : array
        An arbitrarily batched arry of vectors
    vec2 : array
        An arbitrarily batched arry of vectors
    deg : bool (optional)
        Whether angles are returned in degrees or radians. The default is degrees (deg=True).

    Returns
    -------
    angles : array
        A vector of angles with the same leading dimensions of vec1 and vec2.
    """
    v1 = vec1 / np.linalg.norm(vec1, axis=-1)[..., None]
    v2 = vec2 / np.linalg.norm(vec2, axis=-1)[..., None]
    x1 = np.linalg.norm(v1 - v2, axis=-1)
    x2 = np.linalg.norm(v1 + v2, axis=-1)
    alpha = 2.0 * np.arctan2(x1, x2)
    if deg:
        return np.rad2deg(alpha)
    return alpha


def generate_reciprocal_cell(cell, dmin, dtype=np.int32):
    """
    Generate the miller indices of the full P1 reciprocal cell.

    Parameters
    ----------
    cell : tuple, list, np.ndarray of cell parameters, or gemmi.UnitCell
        Unit cell parameters
    dmin : float
        Maximum resolution of the data in Å
    dtype : np.dtype (optional)
        The data type of the returned array. The default is np.int32.


    Returns
    -------
    hkl : np.array(int32)
    """
    hmax, kmax, lmax = cell.get_hkl_limits(dmin)
    hkl = np.meshgrid(
        np.linspace(-hmax, hmax + 1, 2 * hmax + 2, dtype=dtype),
        np.linspace(-kmax, kmax + 1, 2 * kmax + 2, dtype=dtype),
        np.linspace(-lmax, lmax + 1, 2 * lmax + 2, dtype=dtype),
    )
    hkl = np.stack(hkl).reshape((3, -1)).T

    # Remove reflection 0,0,0
    hkl = hkl[np.any(hkl != 0, axis=1)]

    # Remove reflections outside of resolution range
    dHKL = cell.calculate_d_array(hkl).astype("float32")
    hkl = hkl[dHKL >= dmin]

    return hkl


def norm2(array, axis=-1, keepdims=False):
    """Faster version of np.linalg.norm for 3d coords and L2 norm"""
    a2 = np.square(array)
    out = 0.0
    for vec in np.split(a2, a2.shape[axis], axis=axis):
        out += vec
    out = np.sqrt(out)
    if not keepdims:
        out = np.squeeze(out, axis=axis)
    return out


def normalize(array, axis=-1):
    """Normalize a numpy array along a particular axis by dividing by its L2 norm"""
    out = array / norm2(array, axis=axis, keepdims=True)
    return out


class Indexer:
    """
    A class which implements the PinkIndexer algorithm.
    [Gevorkov Y, et al. pinkIndexer – a universal indexer for pink-beam X-ray and electron diffraction snapshots. Acta Cryst A. 2020 Mar 1;76(2):121–31.](https://doi.org/10.1107/S2053273319015559)
    """

    def __init__(self, cell, wav_min, wav_max, float_dtype="float32"):
        """
        Args:
            cell (gemmi.UnitCell): The target cell.
            wav_min (float): Lower wavelength limit in reciprocal Å.
            wav_max (float): Upper wavelength limit in reciprocal Å.
        """
        self.float_dtype = float_dtype
        self.cell = cell
        self.wav_min = wav_min
        self.wav_max = wav_max
        self.B = np.asarray(cell.fractionalization_matrix, dtype=self.float_dtype).T

    def refine_cell_rotation(self, H, scattering_vector, fix="a"):
        """
        Refine rotation and crystal basis from candidate rlp and observed
        scattering vectors.
        Args:
            H (array) : n x 3 array of miller indices.
            scattering_vector (array) : n x 3 array of empirical scattering vectors.
            fix (string) : a comma separated string of cell parameters to fix, possible choices are a,b,c,alpha,beta,gamma
        """
        qhat = normalize(scattering_vector)

        a, b, c, alpha, beta, gamma = self.cell.parameters
        cell_params = np.array([a, b, c, alpha, beta, gamma])
        cell_param_names = ["a", "b", "c", "alpha", "beta", "gamma"]
        fix = fix.split(",")
        fix_mask = np.array(
            [(i in fix) for i in cell_param_names],
        )
        guess = cell_params[~fix_mask]

        def loss_fn(params, return_cell_R_B=False):
            cell = np.copy(cell_params)
            cell[~fix_mask] = params
            cell = gemmi.UnitCell(*cell)

            B = np.array(cell.fractionalization_matrix, dtype=self.float_dtype).T
            rlp = (B @ H.T).T
            rlp_hat = normalize(rlp)

            R, _ = orthogonal_procrustes(rlp_hat, qhat)
            R = R.T
            angles = angle_between((R @ rlp.T).T, qhat)
            loss = angles.mean()

            if return_cell_R_B:
                return loss, cell, R, B
            return loss

        loss_fn(guess)
        result = minimize(loss_fn, guess)
        _, cell, R, B = loss_fn(result.x, True)
        return cell, R

    def index_pink(
        self,
        s0,
        s1,
        max_refls=50,
        rotogram_grid_points=200,
        voxel_grid_points=200,
        fix="a",
    ):
        """
        Args:
            s0 (array): An array of 3 floating point numbers corresponding to the s0 vector. This will be normalized.
            s1 (array): An array of 3 floating point numbers corresponding to the s0 vector. This will be normalized.
            max_refls (int): The maximum number of refls to consider.
            rotogram_grid_points (int): The number of points at which to evaluate the rotograms.
            voxel_grid_points (int): The fineness of the discretization used to quantitate rotograms.
            fix (str): See Indexer.refine_cell_rotation for details.
        """
        s0 = normalize(np.asarray(s0, dtype=self.float_dtype))
        s1 = normalize(np.asarray(s1, dtype=self.float_dtype))

        q = s1 - s0[None, :]
        q_len = norm2(q, keepdims=True)
        qhat = q / q_len

        # Possible resolution range for each observation
        res_min = self.wav_min / q_len.squeeze(-1)
        res_max = self.wav_max / q_len.squeeze(-1)

        if len(res_min) > max_refls:
            dmin = np.sort(res_min)[-max_refls]
        else:
            dmin = res_min.min()

        # Generate the feasible set of reflections from the current geometry
        # These are in cartesian reciprocal space coordinates in the
        # crystal-fixed system
        Hall = generate_reciprocal_cell(self.cell, dmin, dtype=self.float_dtype)
        h = (self.B @ Hall.T).T
        hhat = normalize(h)

        # Remove candidate rlps if incompatible with resolution range of observation
        dall = self.cell.calculate_d_array(Hall)
        mask = (res_max[:, None] >= dall) & (res_min[:, None] <= dall)
        i, j = np.where(mask)

        hhat = hhat[j]
        qhat = qhat[i]

        # mhat bisects hhat and qhat
        m = hhat + qhat
        mhat = normalize(m)

        phi = np.linspace(-2 * np.pi, 0.0, rotogram_grid_points, dtype=self.float_dtype)
        c1 = np.sin(phi / 2.0)
        c2 = -np.cos(phi / 2.0)
        d1 = np.einsum("ad,ad->a", mhat, qhat)
        d2 = np.cross(qhat, mhat)
        theta = 2 * np.arccos(-d1[..., None] * c1)
        num = (
            mhat[..., None, :] * c2[..., :, None] + d2[..., None, :] * c1[..., :, None]
        )
        denom = np.sin(theta[..., None] / 2.0)
        ehat = num / denom
        v = np.arctan(theta / 4)[..., None] * ehat

        # Discretize
        v_grid_size = voxel_grid_points // 2
        rad = np.arctan(np.pi / 4.0)  # Radius of ball in which rotograms live
        _discretized = np.round(v_grid_size * v / rad).astype("int8")
        discretized = _discretized.reshape((-1, 3))

        # This block is optimized for speed. There are lot of nicer, legible ways to
        # do this in numpy but they can be rather slow. Here's one example alternative:
        # idx,counts = np.unique(discretized, axis=0, return_counts=True)
        # max_idx = np.argmax(counts)
        # grid_max = idx[max_idx]
        w = 2 * v_grid_size + 1
        to_id = np.arange(w * w * w, dtype="int32").reshape(
            (w, w, w)
        )  # Converts grid points to a unique integer id
        a, b, c = discretized.T
        flat_ids = to_id[a, b, c]  # Convert 3-D grid values to integers
        counts = np.bincount(flat_ids)
        max_idx = np.argmax(counts)
        brick = np.ones((w, w, w), dtype="int8")
        vrange = np.arange(-v_grid_size, v_grid_size + 1, dtype="int8")
        vrange = np.roll(vrange, v_grid_size + 1)
        grid_max = np.array(
            [
                (vrange[:, None, None] * brick).flatten()[max_idx],
                (vrange[None, :, None] * brick).flatten()[max_idx],
                (vrange[None, None, :] * brick).flatten()[max_idx],
            ]
        )
        # v_max = rad * grid_max / v_grid_size

        # This complicated indexing stuff will figure out the assignment
        # corresponding to the best rotogram voxel.
        flat_assignment = (_discretized == grid_max).all(-1).any(-1)
        Hidx = i[flat_assignment]
        Hobs = Hall[j[flat_assignment]]

        Hout = np.zeros_like(s1)
        Hout[Hidx] = Hobs
        q = s1[Hidx] - s0[None, :]

        # vnorm = norm2(v_max)
        # theta = 4.0 * np.tan(vnorm)

        cell, R = self.refine_cell_rotation(Hobs, q, fix)
        return Hout, cell, R


class PinkIndexer(Strategy):
    """
    A lattice search strategy using the pinkIndexer algorithm.
    For more info, see:
    [Gevorkov Y, et al. pinkIndexer – a universal indexer for pink-beam X-ray and electron diffraction snapshots. Acta Cryst A. 2020 Mar 1;76(2):121–31.](https://doi.org/10.1107/S2053273319015559)
    """

    phil_help = (
        "A lattice search strategy that matches low resolution spots to candidate "
        "indices based on a known unit cell and space group. It supports mono and "
        "polychromatic beams. "
    )

    phil_scope = iotbx.phil.parse(pink_indexer_phil_str)

    def __init__(
        self, target_symmetry_primitive, max_lattices, params=None, *args, **kwargs
    ):
        """Construct PinkIndexer object.
        Args:
            target_symmetry_primitive (cctbx.crystal.symmetry): The target
                crystal symmetry and unit cell
            max_lattices (int): The maximum number of lattice models to find
        """
        super().__init__(params=None, *args, **kwargs)
        self._target_symmetry_primitive = target_symmetry_primitive
        self._max_lattices = max_lattices

        if target_symmetry_primitive is None:
            raise DialsIndexError(
                "Target unit cell and space group must be provided for small_cell"
            )

        target_cell = target_symmetry_primitive.unit_cell()
        if target_cell is None:
            raise ValueError("Please specify known_symmetry.unit_cell")

        self.cell = target_cell
        self.tarsym = target_symmetry_primitive
        self.spacegroup = target_symmetry_primitive.space_group()
        self.wav_min = params.wav_min
        self.wav_max = params.wav_max
        self.max_refls = params.max_refls

    def find_crystal_models(self, reflections, experiments):
        """Find a list of candidate crystal models.
        Args:
            reflections (dials.array_family.flex.reflection_table):
                The found spots centroids and associated data
            experiments (dxtbx.model.experiment_list.ExperimentList):
                The experimental geometry models
        """
        cell = gemmi.UnitCell(*self.cell.parameters())

        wav_min, wav_max = self.wav_min, self.wav_max
        eps = 1e-3
        crystals = []
        for imageset_id, expt in enumerate(experiments):
            beam = expt.beam
            s0 = beam.get_s0()
            wav = beam.get_wavelength()
            if wav_min is None:
                wav_min = wav * (1.0 - eps)
            if wav_max is None:
                wav_max = wav * (1.0 + eps)
            pidxr = Indexer(cell, wav_min, wav_max)
            sel = reflections["imageset_id"] == imageset_id
            refls = reflections.select(sel)
            s1 = np.array(refls["s1"], dtype="float32")
            H, cell, U = pidxr.index_pink(s0, s1, self.max_refls)
            real_a, real_b, real_c = np.array(cell.orthogonalization_matrix).T
            crystal = Crystal(real_a, real_b, real_c, self.spacegroup)
            B = np.array(cell.fractionalization_matrix).T
            UB = U @ B
            A = matrix.sqr(UB.flatten())
            crystal.set_A(A)
            crystals.append(crystal)
        return crystals


#        small_cell_result = small_cell_index_lattice_detail(
#            experiments, reflections, self._params
#        )
#        if not small_cell_result:
#            return []
#
#        ori = small_cell_result[2]
#        direct_matrix = ori.direct_matrix()
#        real_a = direct_matrix[0:3]
#        real_b = direct_matrix[3:6]
#        real_c = direct_matrix[6:9]
#        crystal = Crystal(real_a, real_b, real_c, self._params.small_cell.spacegroup)
#        candidate_crystal_models = [
#            crystal,
#        ]
#
#        # Reset id column to ensure that index assignment later will work properly.
#        # This is a workaround for https://github.com/dials/dials/issues/2485
#        reflections["id"] *= 0
#        reflections["id"] -= 1
#
#        # Additional refinement of the crystal orientation using the spots indexed
#        # by small_cell by superposing ideal relp positions over the observed relps
#        indexed = small_cell_result[3]
#        Bmat = matrix.sqr(crystal.get_B())
#        observed_relps = [(0.0, 0.0, 0.0)]
#        predicted_relps = [(0.0, 0.0, 0.0)]
#        for spot in indexed:
#            predicted_relps.append(Bmat * spot.hkl.ohkl)
#            observed_relps.append(spot.xyz)
#        fit = superpose.least_squares_fit(
#            flex.vec3_double(observed_relps), flex.vec3_double(predicted_relps)
#        )
#        Umat = fit.r
#        crystal.set_A(Umat * Bmat)
#
#        self.candidate_crystal_models = candidate_crystal_models
#        return self.candidate_crystal_models
