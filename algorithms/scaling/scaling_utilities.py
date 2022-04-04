"""
Module of utility functions for scaling.
"""


from __future__ import annotations

import logging
from math import acos

import numpy as np
from scipy.spatial.transform import Rotation

import dxtbx.flumpy as flumpy
from cctbx import miller

from dials.array_family import flex
from dials.util.normalisation import quasi_normalisation as _quasi_normalisation
from dials_scaling_ext import (
    calc_theta_phi,
    create_sph_harm_table,
    rotate_vectors_about_axis,
)

logger = logging.getLogger("dials")

try:
    import platform
    import resource

    def log_memory_usage():
        # getrusage returns kb on linux, bytes on mac
        units_per_mb = 1024
        if platform.system() == "Darwin":
            units_per_mb = 1024 * 1024
        logger.debug(
            "Memory usage: %.1f MB",
            int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / units_per_mb,
        )

except ImportError:

    def log_memory_usage():
        pass


class DialsMergingStatisticsError(Exception):
    """Raised when iotbx merging statistics fails."""

    pass


class BadDatasetForScalingException(Exception):
    """Raised when a selection leaves no further good reflections."""

    pass


class Reasons:
    def __init__(self):
        self.reasons = {}

    def add_reason(self, text, number):
        self.reasons[text] = number

    def __repr__(self):
        reasonlist = [
            f"criterion: {k}, reflections: {v}\n"
            for (k, v) in self.reasons.items()
            if v > 0
        ]
        return "Reflections passing individual criteria:\n" + "".join(reasonlist)


def calc_crystal_frame_vectors(reflection_table, experiment):
    """Calculate the diffraction vectors in the crystal frame."""

    gonio = experiment.goniometer
    fixed_rotation = np.array(gonio.get_fixed_rotation()).reshape(3, 3)
    setting_rotation = np.array(gonio.get_setting_rotation()).reshape(3, 3)
    rotation_axis = np.array(gonio.get_rotation_axis_datum())

    s0c = np.zeros((len(reflection_table), 3))
    s1c = np.zeros((len(reflection_table), 3))
    s0 = np.array(experiment.beam.get_sample_to_source_direction())
    s1 = flumpy.to_numpy(reflection_table["s1"])
    phi = flumpy.to_numpy(
        experiment.scan.get_angle_from_array_index(
            reflection_table["xyzobs.px.value"].parts()[2], deg=False
        )
    )
    # exclude any data that has a bad s1.
    lengths = np.linalg.norm(s1, axis=1)
    non_zero = np.where(lengths > 0.0)
    sel_s1 = s1[non_zero]
    s1n = sel_s1 / lengths[non_zero][:, np.newaxis]
    rotation_matrix = Rotation.from_rotvec(
        phi[non_zero][:, np.newaxis] * rotation_axis
    ).as_matrix()
    R = setting_rotation @ rotation_matrix @ fixed_rotation
    R_inv = np.transpose(R, axes=(0, 2, 1))
    s0c[non_zero] = R_inv @ s0
    # Pairwise matrix multiplication of the arrays of R_inv matrices and s1n vectors
    s1c[non_zero] = np.einsum("ijk,ik->ij", R_inv, s1n)

    reflection_table["s0c"] = flumpy.vec_from_numpy(s0c)
    reflection_table["s1c"] = flumpy.vec_from_numpy(s1c)
    return reflection_table


def align_axis_along_z(alignment_axis, vectors):
    """Rotate the coordinate system such that the exp_rot_axis is along z."""
    if alignment_axis == (0.0, 0.0, 1.0):
        return vectors
    (ux, uy, uz) = alignment_axis
    cross_prod_uz = flex.vec3_double([(uy, -1.0 * ux, 0.0)])
    angle_between_u_z = +1.0 * acos(uz / ((ux**2 + uy**2 + uz**2) ** 0.5))
    phi = flex.double(vectors.size(), angle_between_u_z)
    new_vectors = rotate_vectors_about_axis(cross_prod_uz, vectors, phi)
    return flex.vec3_double(new_vectors)


def sph_harm_table(reflection_table, lmax):
    """Calculate the spherical harmonic table for a spherical
    harmonic absorption correction."""
    theta_phi = calc_theta_phi(reflection_table["s0c"])
    theta_phi_2 = calc_theta_phi(reflection_table["s1c"])
    sph_h_t = create_sph_harm_table(theta_phi, theta_phi_2, lmax)
    return sph_h_t


def quasi_normalisation(reflection_table, experiment):
    """Calculate normalised intensity (Esq) values for reflections, for the purpose
    of selecting subsets based on Esq for scaling. If more involved analyses of
    normalised intensities are needed, then it may be necessary to split this
    procedure to handle acentric and centric reflections separately."""

    logger.info(
        "Calculating normalised intensity values to select a reflection \n"
        "subset for scaling. \n"
    )
    logger.debug(
        "Negative intensities are set to zero for the purpose of \n"
        "calculating mean intensity values for resolution bins. This is to avoid \n"
        "spuriously high E^2 values due to a mean close to zero and should only \n"
        "affect the E^2 values of the highest resolution bins. \n"
    )

    good_refl_sel = ~reflection_table.get_flags(
        reflection_table.flags.bad_for_scaling, all=False
    )
    rt_subset = reflection_table.select(good_refl_sel)

    # Scaling subset is data that has not been flagged as bad or excluded
    miller_set = miller.set(
        crystal_symmetry=experiment.crystal.get_crystal_symmetry(),
        indices=rt_subset["miller_index"],
    )

    # handle negative reflections to minimise effect on mean I values.
    miller_array = miller.array(
        miller_set, data=rt_subset["intensity"], sigmas=rt_subset["variance"] ** 0.5
    )
    if rt_subset.size() <= 10000:
        logger.info(
            """
Insufficient number of reflections (<10000) to calculate normalised intensities.
All reflections will be considered for scaling model determination.
"""
        )
        reflection_table["Esq"] = flex.double(reflection_table.size(), 1.0)
    else:
        normalised_intensities = _quasi_normalisation(miller_array)
        reflection_table["Esq"] = flex.double(reflection_table.size(), 0.0)
        reflection_table["Esq"].set_selected(
            good_refl_sel, normalised_intensities.data()
        )
    return reflection_table


def set_wilson_outliers(reflection_table):
    """Function that takes in a reflection table with 'Esq' and 'centric_flag'
    values and sets an outlier flag depending on a cutoff for p < 1e-6."""

    centric_cutoff = 23.91
    sel1 = reflection_table["centric_flag"]
    sel2 = reflection_table["Esq"] > centric_cutoff  # probability <10^-6
    reflection_table.set_flags(sel1 & sel2, reflection_table.flags.outlier_in_scaling)

    acentric_cutoff = 13.82
    sel1 = ~reflection_table["centric_flag"]
    sel2 = reflection_table["Esq"] > acentric_cutoff  # probability <10^-6
    reflection_table.set_flags(sel1 & sel2, reflection_table.flags.outlier_in_scaling)
    msg = (
        "{0} reflections have been identified as outliers based on their normalised {sep}"
        "intensity values. These are reflections that have a probablity of {sep}"
        "< 10e-6 based on a Wilson distribution (E^2 > {1}, {2} for centric {sep}"
        "and acentric reflections respectively). {sep}"
    ).format(
        reflection_table.get_flags(reflection_table.flags.outlier_in_scaling).count(
            True
        ),
        centric_cutoff,
        acentric_cutoff,
        sep="\n",
    )
    logger.info(msg)
    return reflection_table


def calculate_prescaling_correction(reflection_table):
    """Calculate the multiplicative conversion factor for intensities."""
    conversion = flex.double(reflection_table.size(), 1.0)
    if "lp" in reflection_table:
        conversion *= reflection_table["lp"]
    qe = None
    if "qe" in reflection_table:
        qe = reflection_table["qe"]
    elif "dqe" in reflection_table:
        qe = reflection_table["dqe"]
    if qe:
        inverse_qe = flex.double(reflection_table.size(), 1.0)
        nonzero_qe_sel = qe > 0.0
        good_qe = qe.select(qe > 0.0)
        inverse_qe.set_selected(nonzero_qe_sel.iselection(), 1.0 / good_qe)
        conversion *= inverse_qe
    reflection_table["prescaling_correction"] = conversion
    return reflection_table
