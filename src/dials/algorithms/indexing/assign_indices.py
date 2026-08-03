from __future__ import annotations

import gemmi
import numpy as np
import reciprocalspaceship as rs

from cctbx.array_family import flex
from cctbx.sgtbx import space_group
from cctbx.uctbx import unit_cell

import dials_algorithms_indexing_ext as ext
from dials.algorithms.indexing import DialsIndexError


class AssignIndicesStrategy:
    def __init__(self, d_min=None):
        self._d_min = d_min

    def __call__(self, reciprocal_lattice_vectors):
        raise NotImplementedError()


class AssignIndicesGlobal(AssignIndicesStrategy):
    def __init__(self, tolerance=0.3):
        super().__init__()
        self._tolerance = tolerance

    def __call__(self, reflections, experiments, d_min=None):
        reciprocal_lattice_points = reflections["rlp"]
        reflections["miller_index"] = flex.miller_index(len(reflections), (0, 0, 0))
        if d_min is not None:
            d_spacings = 1 / reciprocal_lattice_points.norms()
            inside_resolution_limit = d_spacings > d_min
        else:
            inside_resolution_limit = flex.bool(reciprocal_lattice_points.size(), True)
        sel = inside_resolution_limit & (reflections["id"] == -1)
        isel = sel.iselection()
        rlps = reciprocal_lattice_points.select(isel)
        refs = reflections.select(isel)
        phi = refs["xyzobs.mm.value"].parts()[2]

        UB_matrices = flex.mat3_double([cm.get_A() for cm in experiments.crystals()])
        imgset_ids = reflections["imageset_id"].select(sel)

        for i_imgset, imgset in enumerate(experiments.imagesets()):
            sel_imgset = imgset_ids == i_imgset

            result = ext.AssignIndices(
                rlps.select(sel_imgset),
                phi.select(sel_imgset),
                UB_matrices,
                tolerance=self._tolerance,
            )

            miller_indices = result.miller_indices()
            crystal_ids = result.crystal_ids()

            expt_ids = flex.int(crystal_ids.size(), -1)
            for i_cryst, cryst in enumerate(experiments.crystals()):
                sel_cryst = crystal_ids == i_cryst
                for i_expt in experiments.where(crystal=cryst, imageset=imgset):
                    expt_ids.set_selected(sel_cryst, i_expt)
                    if experiments[i_expt].identifier:
                        reflections.experiment_identifiers()[i_expt] = experiments[
                            i_expt
                        ].identifier

            reflections["miller_index"].set_selected(
                isel.select(sel_imgset), miller_indices
            )
            reflections["id"].set_selected(isel.select(sel_imgset), expt_ids)
            reflections.set_flags(
                reflections["miller_index"] != (0, 0, 0), reflections.flags.indexed
            )
            reflections["id"].set_selected(reflections["miller_index"] == (0, 0, 0), -1)


class AssignIndicesLocal(AssignIndicesStrategy):
    def __init__(
        self, d_min=None, epsilon=0.05, delta=8, l_min=0.8, nearest_neighbours=20
    ):
        super().__init__()
        self._epsilon = epsilon
        self._delta = delta
        self._l_min = l_min
        self._nearest_neighbours = nearest_neighbours

    def __call__(self, reflections, experiments, d_min=None):
        from libtbx.math_utils import nearest_integer as nint
        from scitbx import matrix

        reciprocal_lattice_points = reflections["rlp"]
        if "miller_index" not in reflections:
            reflections["miller_index"] = flex.miller_index(len(reflections))
        if d_min is not None:
            d_spacings = 1 / reciprocal_lattice_points.norms()
            inside_resolution_limit = d_spacings > d_min
        else:
            inside_resolution_limit = flex.bool(reciprocal_lattice_points.size(), True)
        sel = inside_resolution_limit & (reflections["id"] == -1)
        isel = sel.iselection()
        rlps = reciprocal_lattice_points.select(isel)
        refs = reflections.select(isel)
        phi = refs["xyzobs.mm.value"].parts()[2]

        if len(rlps) <= self._nearest_neighbours:
            raise DialsIndexError(
                "index_assignment.local.nearest_neighbour must be smaller than the number of accepted reflections (%d)"
                % len(rlps)
            )

        UB_matrices = flex.mat3_double([cm.get_A() for cm in experiments.crystals()])
        imgset_ids = reflections["imageset_id"].select(sel)

        for i_imgset, imgset in enumerate(experiments.imagesets()):
            sel_imgset = imgset_ids == i_imgset

            result = ext.AssignIndicesLocal(
                rlps.select(sel_imgset),
                phi.select(sel_imgset),
                UB_matrices,
                epsilon=self._epsilon,
                delta=self._delta,
                l_min=self._l_min,
                nearest_neighbours=self._nearest_neighbours,
            )
            miller_indices = result.miller_indices()
            crystal_ids = result.crystal_ids()
            hkl = miller_indices.as_vec3_double().iround()

            assert miller_indices.select(crystal_ids < 0).all_eq((0, 0, 0))

            expt_ids = flex.int(crystal_ids.size(), -1)
            for i_cryst, cryst in enumerate(experiments.crystals()):
                sel_cryst = crystal_ids == i_cryst
                for i_expt in experiments.where(crystal=cryst, imageset=imgset):
                    expt_ids.set_selected(sel_cryst, i_expt)
                    if experiments[i_expt].identifier:
                        reflections.experiment_identifiers()[i_expt] = experiments[
                            i_expt
                        ].identifier

            for i_cryst in set(crystal_ids):
                if i_cryst < 0:
                    continue

                A = matrix.sqr(experiments[i_cryst].crystal.get_A())
                A_inv = A.inverse()

                cryst_sel = crystal_ids == i_cryst
                rlp_sel = rlps.select(cryst_sel)
                hkl_sel = hkl.select(cryst_sel).as_vec3_double()

                d_sel = 1 / rlp_sel.norms()
                d_perm = flex.sort_permutation(d_sel, reverse=True)

                hf_0 = A_inv * rlp_sel[d_perm[0]]
                h_0 = matrix.col([nint(j) for j in hf_0.elems])
                offset = h_0 - matrix.col(hkl_sel[d_perm[0]])
                # print "offset:", offset.elems

                h = hkl_sel + flex.vec3_double(hkl_sel.size(), offset.elems)

                refs["miller_index"].set_selected(
                    cryst_sel, flex.miller_index(list(h.iround()))
                )
                refs["id"].set_selected(cryst_sel, i_cryst)

            crystal_ids.set_selected(crystal_ids < 0, -1)
            refs["id"] = crystal_ids
            refs["miller_index"].set_selected(crystal_ids < 0, (0, 0, 0))

            reflections["miller_index"].set_selected(isel, refs["miller_index"])
            reflections["id"].set_selected(isel, refs["id"])
            reflections.set_flags(
                reflections["miller_index"] != (0, 0, 0), reflections.flags.indexed
            )


class AssignIndicesPolychromatic(AssignIndicesStrategy):
    """
    An object to assign Miller indices to a polychromatic Laue still image.
    """

    def __init__(self, lam_min, lam_max, lam_peak, d_min, spacegroup="1"):
        """
        Parameters:
            s0 (np.ndarray): A 3-vector indicating the direction of the incoming beam wavevector.
            s1 (np.ndarray): n x 3 array indicating the direction of the scattered beam wavevector.
            cell (Union[Tuple, gemmi.UnitCell]): A tuple or list of unit cell params (a, b, c, alpha, beta, gamma) or a gemmi.UnitCell object.
            R (np.ndarray): A 3x3 rotation matrix corresponding to the crystal orientation for the frame.
            lam_min (float): The lower end of the wavelength range of the beam.
            lam_max (float): The upper end of the wavelength range of the beam.
            lam_peak (float): The wavelength with peak intensity.
            d_min (float): The maximum resolution of the model.
            spacegroup (str): Anything that the gemmi.SpaceGroup constructor understands.
        """
        super().__init__()

        self.lam_min = lam_min
        self.lam_max = lam_max
        self.lam_peak = lam_peak
        self.d_min = d_min

    def __call__(self, reflections, experiments, d_min=None):
        """
        Runs indexing, wavelength assignment, and orientation refinement jointly.
        """
        # Get experiment variables
        s0 = np.array(experiments[0].beam.get_s0())
        s1 = reflections["s1"].as_numpy_array()
        cryst = experiments[0].crystal
        self.R = np.asarray(cryst.get_U()).reshape(3, 3)
        cell_params = cryst.get_unit_cell().parameters()
        self.cell = gemmi.UnitCell(*cell_params)
        self.spacegroup = gemmi.SpaceGroup(
            cryst.get_space_group().type().universal_hermann_mauguin_symbol()
        )

        self.s0 = s0 / np.linalg.norm(s0)
        self.B = np.array(self.cell.frac.mat).T

        # Initialize the full reciprocal grid
        hmax, kmax, lmax = self.cell.get_hkl_limits(self.d_min)
        Hall = (
            np.mgrid[
                -hmax : hmax + 1 : 1.0,
                -kmax : kmax + 1 : 1.0,
                -lmax : lmax + 1 : 1.0,
            ]
            .reshape((3, -1))
            .T
        )
        Hall = Hall[np.any(Hall != 0, axis=1)]
        d = self.cell.calculate_d_array(Hall)
        Hall = Hall[d >= self.d_min]

        # Remove any systematic absences in the space group
        Hall = Hall[~rs.utils.is_absent(Hall, self.spacegroup)]
        self.Hall = Hall

        # Initialize class variables needed for assignment
        self._s1 = s1 / np.linalg.norm(s1, axis=-1)[:, None]
        self._qobs = self._s1 - self.s0
        self._qpred = np.zeros_like(self._s1)
        self._H = np.zeros_like(self._s1)
        self._wav = np.zeros(len(self._H))
        self._harmonics = np.zeros(len(self._s1), dtype=bool)
        self._inliers = np.ones(len(self._s1), dtype=bool)

        self.assign()
        for j in range(5):
            self.reset_inliers()
            self.update_rotation()
            self.assign()
            self.reject_outliers()
            self.update_rotation()
            self.assign()

        # Update s1 based on new wavelengths
        np.divide(
            self._s1, self._wav[:, None], out=self._s1, where=self._wav[:, None] != 0
        )

        # Reset crystal parameters based on new geometry
        experiments[0].crystal.set_U(self.R.flatten())
        experiments[0].crystal.set_A(self.RB.flatten())
        experiments[0].crystal.set_B(self.B.flatten())
        experiments[0].crystal.set_space_group(space_group(self.spacegroup.hall))
        experiments[0].crystal.set_unit_cell(unit_cell(self.cell.parameters))

        # Get wavelengths
        spot_wavelengths = np.asarray(self._wav.tolist())

        # Write data to reflections
        reflections["s1"] = flex.vec3_double(self._s1)
        reflections["miller_index"] = flex.miller_index(self._H.astype("int").tolist())
        reflections["harmonics"] = flex.bool(self._harmonics.tolist())
        reflections["wavelength"] = flex.double(spot_wavelengths)
        reflections["id"] = flex.int([0] * len(reflections))

        # Filter reflections by spectrum
        all_wavelengths = reflections["wavelength"].as_numpy_array()
        keep = np.logical_and(
            all_wavelengths >= self.lam_min,
            all_wavelengths <= self.lam_max,
        )
        reflections["id"].set_selected(flex.bool(~keep), -1)

        # Set flags for indexed reflections
        reflections.set_flags(
            reflections["miller_index"] != (0, 0, 0), reflections.flags.indexed
        )
        reflections["id"].set_selected(reflections["miller_index"] == (0, 0, 0), -1)

    @property
    def RB(self):
        """
        Calculates the product RB.

        Returns:
            np.ndarray: The product of the rotation matrix (R) and the fractionalization matrix (B).
        """
        return self.R @ self.B

    @property
    def s1(self):
        """
        np.ndarray: The normalized direction of the scattered beam wavevector for inlying reflections.
        """
        return self._s1[self._inliers]

    @property
    def qobs(self):
        """
        np.ndarray: The observed q vectors for inlying reflections.
        """
        return self._qobs[self._inliers]

    @property
    def qpred(self):
        """
        np.ndarray: The predicted q vectors for inlying reflections.
        """
        return self._qpred[self._inliers]

    @property
    def H(self):
        """
        np.ndarray: The Miller indices for inlying reflections.
        """
        return self._H[self._inliers]

    @property
    def wav(self):
        """
        np.ndarray: The wavelengths associated with inlying reflections.
        """
        return self._wav[self._inliers]

    @property
    def harmonics(self):
        """
        np.ndarray: Boolean array indicating whether inlying reflections correspond to harmonic reflections.
        """
        return self._harmonics[self._inliers]

    # <-- setters that operate on the currently inlying set
    def set_qpred(self, qpred):
        """
        Set the predicted q vectors for inlying reflections.

        Parameters:
            qpred (np.ndarray): Predicted q vectors.
        """
        self._qpred[self._inliers] = qpred

    def set_H(self, H):
        """
        Set the Miller indices for inlying reflections.

        Parameters:
            H (np.ndarray): Miller indices.
        """
        self._H[self._inliers] = H
        self._H[~self._inliers] = 0.0

    def set_wav(self, wav):
        """
        Set the wavelengths for inlying reflections.

        Parameters:
            wav (np.ndarray): Wavelengths.
        """
        self._wav[self._inliers] = wav
        self._wav[~self._inliers] = 0.0

    def set_inliers(self, inliers):
        """
        Set the inliers for the current set of reflections.

        Parameters:
            inliers (np.ndarray): Boolean array indicating inliers.
        """
        self._inliers[self._inliers] = inliers

    def set_harmonics(self, harmonics):
        """
        Set whether inlying reflections correspond to harmonics.

        Parameters:
            harmonics (np.ndarray): Boolean array indicating harmonics.
        """
        self._harmonics[self._inliers] = harmonics

    # --> setters that operate on the currently inlying set

    def reset_inliers(self):
        """Reset all reflections as inliers."""
        self._inliers = np.ones(len(self._inliers), dtype=bool)

    def reject_outliers(self, nstd=10.0):
        """
        Update the list of inliers based on robust statistical measures.

        Parameters:
            nstd (float): Number of standard deviations from the median to consider as inliers.

        This method uses Minimum Covariance Determinant (MCD) to robustly estimate the covariance matrix
        of the concatenated observed and predicted q vectors. It then considers inliers as the points
        within a specified number of standard deviations from the median Mahalanobis distance.

        The list of inliers is updated accordingly.
        """
        from sklearn.covariance import MinCovDet

        X = np.concatenate((self.qobs, self.qpred * self.wav[:, None]), axis=-1)
        dist = MinCovDet().fit(X).dist_
        self.set_inliers(dist <= nstd**2.0)

    def assign(self):
        """
        Assign miller indices to the inlier reflections.

        This method updates the following attributes:
        - self.H: Miller indices associated with the assigned reflections.
        - self.wav: Wavelengths associated with the assigned reflections.
        - self.qpred: Predicted scattering vectors associated with the assigned reflections.
        - self.harmonics: Boolean array indicating whether the assigned reflections are harmonics.

        The assignment is performed by solving the linear sum assignment problem using scipy.optimize.linear_sum_assignment.
        The cost matrix is computed based on the angular distance between observed and predicted scattering vectors.

        The feasible set of reflections is determined from the reciprocal lattice points within the specified geometry.
        Reflections with duplicated scattering vectors are removed, and then the assignment is performed.

        This method is essential for updating the information of assigned reflections, allowing subsequent steps in
        Laue indexing procedures to utilize the assigned miller indices and associated parameters.
        """
        # Generate the feasible set of reflections from the current geometry
        Hall = self.Hall
        qall = (self.RB @ Hall.T).T
        feasible = (
            np.linalg.norm(qall + self.s0 / self.lam_min, axis=-1) < 1 / self.lam_min
        ) & (np.linalg.norm(qall + self.s0 / self.lam_max, axis=-1) > 1 / self.lam_max)
        Hall = Hall[feasible]
        qall = qall[feasible]

        # Keep track of harmonics in the feasible set
        Raypred = hkl2ray(Hall)
        _, idx, counts = np.unique(
            Raypred, return_index=True, return_counts=True, axis=0
        )
        unique_rays, counts = np.unique(Raypred, return_counts=True, axis=0)
        harmonic_rays = counts > 1

        # Keep only harmonics closest to peak wavelength in feasible set
        harmonics = np.zeros(len(Hall), dtype=bool)
        to_keep = np.ones(len(Hall), dtype=bool)
        for i, ray in enumerate(unique_rays[harmonic_rays]):
            idh = np.where((Raypred == ray).all(axis=1))[0]
            if len(idh) == 1:  # Not a harmonic
                to_keep[idh] = True
            else:
                harmonics[idh] = True
                qharm = qall[idh]
                harmonic_wavelengths = (
                    -2.0 * (self.s0 * qharm).sum(-1) / (qharm * qharm).sum(-1)
                )
                kept_harmonic = np.argmin(np.abs(harmonic_wavelengths - self.lam_peak))
                to_keep[idh] = False
                to_keep[idh[kept_harmonic]] = True

        # Remove non-optimal harmonics from the feasible set
        Hall = Hall[to_keep]
        qall = qall[to_keep]
        harmonics = harmonics[to_keep]

        dmat = rs.utils.angle_between(self.qobs[..., None, :], qall[None, ..., :])
        cost = dmat

        from scipy.optimize import linear_sum_assignment

        ido, idx = linear_sum_assignment(cost)

        # Update appropriate variables
        H = Hall[idx]
        qpred = qall[idx]
        harmonics = harmonics[idx]

        # Set all attributes to match the current assignment
        self.set_H(H)
        self.set_qpred(qpred)
        self.set_harmonics(harmonics)

        # wav_pred = -2.*(self.s0 * qpred).sum(-1) / (qpred*qpred).sum(-1)
        with np.errstate(divide="ignore"):
            wav_obs = np.linalg.norm(self.qobs, axis=-1) / np.linalg.norm(
                self.qpred, axis=-1
            )
        self.set_wav(wav_obs)

    def update_rotation(self):
        """Update the rotation matrix (self.R) based on the inlying reflections."""
        from scipy.linalg import orthogonal_procrustes

        misset, _ = orthogonal_procrustes(
            self.qobs,
            self.qpred * self.wav[:, None],
        )
        self.R = misset @ self.R


def hkl2ray(hkl, wavelength=None):
    """
    Convert a miller index to the shortest member of its central ray.
    Optionally, adjust its wavelength accordingly.

    Parameters:
        hkl (np.array): `n x 3` array of miller indices. The dtype must be interpretable as an integer.
        wavelength (Optional[np.array]): Length `n` array of wavelengths corresponding to each miller index.

    Returns:
        reduced_hkl (np.array): The miller index of the shortest vector on the same central ray as the original hkl.
        reduced_wavelength (Optional[np.array]): The wavelengths corresponding to reduced_hkl.
    """
    gcd = np.gcd.reduce(hkl.astype(int), axis=-1)
    if wavelength is not None:
        return hkl / gcd[..., None], wavelength * gcd
    else:
        return hkl / gcd[..., None]
