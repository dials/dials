from __future__ import annotations

import logging

import gemmi
import numpy as np
from scipy.spatial.transform.rotation import Rotation

import iotbx.phil
from dxtbx.model import Crystal

from dials.algorithms.indexing import DialsIndexError

from .strategy import Strategy

logger = logging.getLogger(__name__)

pink_indexer_phil_str = """
pink_indexer
    .expert_level = 1
{
    max_refls = 50
        .type = int(value_min=10)
        .help = "Maximum number of reflections to consider indexing"
    wavelength = None
        .type = float(value_min=0.)
        .help = "The peak wavelength"
    percent_bandwidth = 1.
        .type = float(value_min=0.)
        .help = "The percent bandwidth used to calculate the wavelength range for indexing. The wavelength range is defined (wavelength - wavelength*percent_bandwidth/200, wavelength + wavelength*percent_bandwidth/200). This parameter also reflects the uncertainty of the supplied cell constants with larger values appropriate for less certain unit cells."
    rotogram_grid_points = 180
        .type = int(value_min=10, value_max=1000)
        .help = "Number of points at which to evaluate the angle search for each rlp-observation pair"
    voxel_grid_points=150
        .type = int(value_min=10, value_max=1000)
        .help = "Controls the number of voxels onto which the rotograms are discretized"
    min_lattices=1
        .type = int(value_min=1, value_max=100)
        .help = "The minimum number of candidate lattices to generate."
}
"""


def rotvec_to_quaternion(rotvec, deg=False, eps=1e-32):
    """
    Convert rotation vector(s) to quaternion(s).
    """
    alpha = norm2(rotvec)
    ax = rotvec / np.where(alpha == 0.0, np.inf, alpha)[..., None]
    if deg:
        alpha = np.deg2rad(alpha)
    a2 = 0.5 * alpha
    w = np.cos(a2)
    x = np.sin(a2) * ax[..., 0]
    y = np.sin(a2) * ax[..., 1]
    z = np.sin(a2) * ax[..., 2]
    return np.stack((w, x, y, z), axis=-1)


def quaternion_multiply(a, b):
    """
    Multiply two quaternions, return a*b
    """
    # Expand a
    aw = a[..., 0]
    ax = a[..., 1]
    ay = a[..., 2]
    az = a[..., 3]

    # Expand b
    bw = b[..., 0]
    bx = b[..., 1]
    by = b[..., 2]
    bz = b[..., 3]

    # Calculate result
    w = aw * bw - ax * bx - ay * by - az * bz
    x = aw * bx + ax * bw + ay * bz - az * by
    y = aw * by - ax * bz + ay * bw + az * bx
    z = aw * bz + ax * by - ay * bx + az * bw

    return np.dstack((w, x, y, z))


def norm2(array, axis=-1, keepdims=False):
    """Faster version of np.linalg.norm for the L2 norm in lower dimensions."""
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
    v1 = normalize(vec1, axis=-1)
    v2 = normalize(vec2, axis=-1)
    x1 = norm2(v1 - v2, axis=-1)
    x2 = norm2(v1 + v2, axis=-1)
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


class Indexer:
    """
    A class which implements the PinkIndexer algorithm.
    [Gevorkov Y, et al. pinkIndexer – a universal indexer for pink-beam X-ray and electron diffraction snapshots. Acta Cryst A. 2020 Mar 1;76(2):121–31.](https://doi.org/10.1107/S2053273319015559)
    """

    def __init__(self, cell, wavelength, bandwidth, float_dtype="float32"):
        """
        Args:
            cell (gemmi.UnitCell): The target cell.
            wavelength (float): Peak wavelength in Å.
            bandwidth (float): The percentage bandwidth used to calculate a wavelength range
        """
        self.float_dtype = float_dtype
        self.cell = cell
        self.wav_peak = wavelength
        self.wav_min = wavelength - bandwidth * wavelength / 200.0
        self.wav_max = wavelength + bandwidth * wavelength / 200.0
        self.B = np.asarray(cell.fractionalization_matrix, dtype=self.float_dtype).T

    def index_pink(
        self,
        s0,
        s1,
        max_refls=50,
        rotogram_grid_points=200,
        voxel_grid_points=200,
        int_dtype="uint8",
        float_dtype="float32",
        min_lattices=1,
        dilate_r=None,
    ):
        """
        Args:
            s0 (array): An array of 3 floating point numbers corresponding to the s0 vector. This will be normalized.
            s1 (array): An n x 3 array floating point numbers corresponding to the s1 vectors. This will be normalized.
            max_refls (int, optional): The maximum number of refls to consider.
            rotogram_grid_points (int, optional): The number of points at which to evaluate the rotograms.
            voxel_grid_points (int, optional): The fineness of the discretization used to quantitate rotograms.
            int_dtype (str or dtype, optional): the dtype to use for the voxel grid
            float_dtype (str or dtype, optional): the dtype to use for floating point values
            min_lattices (int, optional): the minimum number of candidate lattices returned by this function
            dilate_r (float, optional): optionally dilate the voxel grid by a kernel with this radius in pixels.

        Returns:
            iter_UB: an interable of crystal bases as numpy arrays
        """
        s0_hat = normalize(np.asarray(s0, dtype=self.float_dtype))
        s1_hat = normalize(np.asarray(s1, dtype=self.float_dtype))

        q = s1_hat - s0_hat[None, :]
        q_len = norm2(q, keepdims=True)  # length of q if lambda is 1A

        # Possible resolution range for each observation considering wavelength range
        res_min = self.wav_min / q_len.squeeze(-1)
        res_max = self.wav_max / q_len.squeeze(-1)

        in_range = None
        # Truncate the resolution range if there are too many reflections
        if len(res_min) > max_refls:
            dmin = np.sort(res_min)[-max_refls]
            in_range = res_min >= dmin
            s1_hat = s1_hat[in_range]
            s1 = s1[in_range]
            q = q[in_range]
            q_len = q_len[in_range]
            res_min = res_min[in_range]
            res_max = res_max[in_range]
        else:
            dmin = res_min.min()

        # Normalized q vector
        qhat = normalize(q)

        # Generate the feasible set of reflections from the current geometry
        # These are in cartesian reciprocal space coordinates in the
        # crystal-fixed system
        Hall = generate_reciprocal_cell(self.cell, dmin, dtype=self.float_dtype)
        h = (self.B @ Hall.T).T
        hhat = normalize(h)

        # Remove candidate rlps if incompatible with resolution range of observation
        dall = self.cell.calculate_d_array(Hall)
        mask = (dall <= res_max[:, None]) & (dall >= res_min[:, None])
        i, j = np.where(mask)

        hhat = hhat[j]
        qhat = qhat[i]

        # mhat bisects hhat and qhat
        m = hhat + qhat
        mhat = normalize(m)

        # construct a rotation that maps hhat onto qhat
        Rm = rotvec_to_quaternion(np.pi * mhat)

        # construct rotations about qhat
        phimin = -np.pi
        phimax = np.pi
        phi = np.linspace(
            phimin, phimax, rotogram_grid_points + 1, dtype=self.float_dtype
        )[:-1]
        rotvec = qhat[:, None, :] * phi[None, :, None]
        Rq = rotvec_to_quaternion(rotvec)

        # combine mhat and qhat rotations into general rotation
        quat = quaternion_multiply(Rq, Rm[:, None, :])
        flat_quat = quat.reshape((-1, 4))
        rotvec = (
            Rotation.from_quat(flat_quat).as_rotvec().reshape(quat.shape[:-1] + (-1,))
        )
        theta = norm2(rotvec)
        axis = rotvec / theta[..., None]

        # Scaling the general rotvec norm by this factor makes the discretization more uniform
        # without this, the rotations would be biased toward higher angles
        scales = np.arctan(theta / 4.0)
        scaled_rotvec = axis * scales[..., None]

        # Discretize scaled rotvecs
        scale_max = np.arctan(np.pi / 4.0)
        bins = np.linspace(-scale_max, scale_max, voxel_grid_points)
        # This is how you would calculate the bin centers if you cared to extract the rotation matrix directly
        # bin_centers = np.concatenate(
        #    (bins[[0]], 0.5 * (bins[1:] + bins[:-1]), bins[[-1]])
        # )
        idx = np.digitize(scaled_rotvec, bins)

        # Map discretized rotvecs into voxel grid
        n = voxel_grid_points + 1
        voxels = np.zeros((n, n, n), dtype=int_dtype)
        np.add.at(voxels, tuple(idx.transpose((2, 0, 1))), 1)

        # Optionally dilate the voxel grid
        if dilate_r is not None:
            # PinkIndexer does a little dilation to help avoid overfitting
            # I'm implementing this using a radially symmetric convolution
            # In the original paper, they use a cubic kernel ones((3, 3, 3))
            voxels = voxels.astype(float_dtype)
            x = np.arange(voxels.shape[0], dtype=float_dtype)
            x = np.square(x - x.mean())
            kernel = np.exp(
                -(x[:, None, None] + x[None, :, None] + x[None, None, :])
                / np.square(dilate_r)
            )
            from scipy.signal import fftconvolve

            voxels = fftconvolve(voxels, kernel, mode="same")

        # Possible solutions are voxels with the highest density
        cutoff = np.sort(voxels.flatten())[-min_lattices]
        peaks = np.column_stack(np.where(voxels >= cutoff))
        for peak in peaks:
            assignment = np.zeros_like(mask)
            assignment[i, j] = (idx == peak).all(-1).any(-1)
            refl_id, miller_id = np.where(assignment)

            # Use a bootstrap approach to estimate the UB matrix
            n = len(refl_id)
            num_bootstraps = 100
            bs_idx = np.random.choice(n, (num_bootstraps, n))

            refl_id = refl_id[bs_idx]
            miller_id = miller_id[bs_idx]

            h = Hall[miller_id]
            s = s1_hat[refl_id]

            # Assign everything to be the nominal wavelength
            wav = self.wav_peak
            k = np.reciprocal(wav)
            UB = (
                k
                * (s - s0_hat).transpose(0, 2, 1)
                @ np.linalg.pinv(h.transpose(0, 2, 1))
            )
            yield UB.mean(0)


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
        """
        Args:
            target_symmetry_primitive (cctbx.crystal.symmetry): The target crystal symmetry and unit cell
            max_lattices (int): The maximum number of lattice models to find
            params (phil,optional): Phil params
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
        self.wavelength = params.wavelength
        self.percent_bandwidth = params.percent_bandwidth
        self.max_refls = params.max_refls
        self.rotogram_grid_points = params.rotogram_grid_points
        self.voxel_grid_points = params.voxel_grid_points
        self.min_lattices = params.min_lattices

    def find_crystal_models(self, reflections, experiments):
        """Find a list of candidate crystal models.
        Args:
            reflections (dials.array_family.flex.reflection_table):
                The found spots centroids and associated data
            experiments (dxtbx.model.experiment_list.ExperimentList):
                The experimental geometry models
        """
        # This is a workaround for https://github.com/dials/dials/issues/2485
        reflections["id"] *= 0
        reflections["id"] -= 1

        if len(experiments) < 1:
            msg = "pink_indexer received an experimentlist with length > 1. To use the pink_indexer method, you must set joint_indexing=False when you call dials.index"
            raise ValueError(msg)

        expt = experiments[0]
        refls = reflections

        cell = gemmi.UnitCell(*self.cell.parameters())

        beam = expt.beam
        s0 = beam.get_s0()
        wav, bw = self.wavelength, self.percent_bandwidth
        if wav is None:
            wav = beam.get_wavelength()
        if bw is None:
            bw = 1.0
        pidxr = Indexer(cell, wav, bw)
        s1 = np.array(refls["s1"], dtype="float32")
        self.candidate_crystal_models = []
        for UB in pidxr.index_pink(
            s0,
            s1,
            self.max_refls,
            rotogram_grid_points=self.rotogram_grid_points,
            voxel_grid_points=self.voxel_grid_points,
            min_lattices=self.min_lattices,
        ):
            real_a, real_b, real_c = np.linalg.inv(UB.astype("double"))
            crystal = Crystal(real_a, real_b, real_c, self.spacegroup)
            self.candidate_crystal_models.append(crystal)
        return self.candidate_crystal_models
