"""
Module of utility functions for scaling.
"""


import logging
from math import acos

import scitbx
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
    fixed_rotation = scitbx.matrix.sqr(gonio.get_fixed_rotation())
    setting_rotation = scitbx.matrix.sqr(gonio.get_setting_rotation())
    rotation_axis = scitbx.matrix.col(gonio.get_rotation_axis_datum())

    s0c = flex.vec3_double(reflection_table.size(), (0, 0, 0))
    s1c = flex.vec3_double(reflection_table.size(), (0, 0, 0))
    # we want sample to source direction.
    s0 = tuple(-1.0 * i for i in experiment.beam.get_unit_s0())
    # exclude any data that has a bad s1.
    lengths = flex.double(scitbx.matrix.col(v).length() for v in reflection_table["s1"])
    sel = lengths > 0.0
    sel_s1 = reflection_table["s1"].select(sel)
    sel_z = reflection_table["xyzobs.px.value"].parts()[2].select(sel)
    s1n = sel_s1.each_normalize()
    for z, s1, index in zip(sel_z, s1n, sel.iselection()):
        phi = experiment.scan.get_angle_from_array_index(z, deg=True)
        rotation_matrix = rotation_axis.axis_and_angle_as_r3_rotation_matrix(
            phi, deg=True
        )
        R_inv = (setting_rotation * rotation_matrix * fixed_rotation).inverse()
        s0c[index] = R_inv * s0
        s1c[index] = R_inv * s1

    reflection_table["s0c"] = s0c
    reflection_table["s1c"] = s1c

    return reflection_table


def align_axis_along_z(alignment_axis, vectors):
    """Rotate the coordinate system such that the exp_rot_axis is along z."""
    if alignment_axis == (0.0, 0.0, 1.0):
        return vectors
    (ux, uy, uz) = alignment_axis
    cross_prod_uz = flex.vec3_double([(uy, -1.0 * ux, 0.0)])
    angle_between_u_z = +1.0 * acos(uz / ((ux ** 2 + uy ** 2 + uz ** 2) ** 0.5))
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
