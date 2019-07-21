"""
Module of utility functions for scaling.
"""

from __future__ import absolute_import, division, print_function

import copy
import logging
from math import pi, acos
from time import time

from dials.array_family import flex
from cctbx import miller
from cctbx import uctbx
from dxtbx.model.experiment_list import ExperimentListDumper
from dials_scaling_ext import (
    create_sph_harm_table,
    calc_theta_phi,
    rotate_vectors_about_axis,
)


logger = logging.getLogger("dials")

try:
    import resource
    import platform

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


class Reasons(object):
    def __init__(self):
        self.reasons = {}

    def add_reason(self, text, number):
        self.reasons[text] = number

    def __repr__(self):
        reasonlist = [
            "criterion: %s, reflections: %s\n" % (k, v)
            for (k, v) in self.reasons.items()
            if v > 0
        ]
        return "Reflections passing individual criteria:\n" + "".join(reasonlist)


def save_experiments(experiments, filename):
    """Save the experiments json."""
    st = time()
    logger.info("Saving the experiments to %s", filename)
    dump = ExperimentListDumper(experiments)
    with open(filename, "w") as outfile:
        outfile.write(dump.as_json(split=True))
    logger.info("Time taken: %g", (time() - st))


def save_reflections(reflection_table, filename):
    """Save the scaled reflections."""
    st = time()
    to_del = [
        "variance",
        "intensity",
        "s0",
        "s0c",
        "s1c",
        "prescaling_correction",
        "batch",
    ]
    for col in to_del:
        try:
            del reflection_table[col]
        except KeyError:
            pass
    logger.info("Saving the scaled reflections to %s", filename)
    reflection_table.as_pickle(filename)
    logger.info("Time taken: %g", (time() - st))


"""def calc_sigmasq(jacobian_transpose, var_cov):
  sigmasq = flex.float([])
  #note: must be a faster way to do this next bit? - in c++?
  for col in jacobian_transpose.cols(): #iterating over reflections
    a = flex.double(col.as_dense_vector())
    var = (a * var_cov) * a
    sigmasq.append(flex.sum(var))
  return sigmasq.as_double()"""


def calc_crystal_frame_vectors(reflection_table, experiments):
    """Calculate the diffraction vectors in the crystal frame."""
    reflection_table["s0"] = flex.vec3_double(
        [experiments.beam.get_s0()] * len(reflection_table)
    )
    rot_axis = flex.vec3_double([experiments.goniometer.get_rotation_axis()])
    angles = reflection_table["phi"] * -1.0 * pi / 180  # want to do an inverse rot.
    reflection_table["s1c"] = rotate_vectors_about_axis(
        rot_axis, reflection_table["s1"], angles
    )
    reflection_table["s0c"] = rotate_vectors_about_axis(
        rot_axis, reflection_table["s0"], angles
    )
    reflection_table["s1c"] = align_rotation_axis_along_z(
        rot_axis, reflection_table["s1c"]
    )
    reflection_table["s0c"] = align_rotation_axis_along_z(
        rot_axis, reflection_table["s0c"]
    )
    return reflection_table


def align_rotation_axis_along_z(exp_rot_axis, vectors):
    """Rotate the coordinate system such that the exp_rot_axis is along z."""
    if list(exp_rot_axis) == [(0.0, 0.0, 1.0)]:
        return vectors
    (ux, uy, uz) = exp_rot_axis[0][0], exp_rot_axis[0][1], exp_rot_axis[0][2]
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
    rt_subset["intensity_for_norm"] = copy.deepcopy(rt_subset["intensity"])
    rt_subset["intensity_for_norm"].set_selected(rt_subset["intensity"] < 0.0, 0.0)
    miller_array = miller.array(miller_set, data=rt_subset["intensity_for_norm"])
    n_refl = rt_subset.size()

    # set up binning objects
    if n_refl > 20000:
        n_refl_shells = 20
    elif n_refl > 15000:
        n_refl_shells = 15
    elif n_refl > 10000:
        n_refl_shells = 10
    else:
        logger.info(
            "No normalised intensity values were calculated, as an insufficient\n"
            "number of reflections were detected. All normalised intensity \n"
            "values will be set to 1 to allow use in scaling model determination. \n"
        )
        reflection_table["Esq"] = flex.double(reflection_table.size(), 1.0)
        return reflection_table

    d_star_sq = miller_array.d_star_sq().data()
    d_star_sq_min = flex.min(d_star_sq)
    d_star_sq_max = flex.max(d_star_sq)
    span = d_star_sq_max - d_star_sq_min

    relative_tolerance = 1e-6
    d_star_sq_max += span * relative_tolerance
    d_star_sq_min -= span * relative_tolerance
    step = (d_star_sq_max - d_star_sq_min) / n_refl_shells
    miller_array.setup_binner_d_star_sq_step(
        auto_binning=False,
        d_max=uctbx.d_star_sq_as_d(d_star_sq_max),
        d_min=uctbx.d_star_sq_as_d(d_star_sq_min),
        d_star_sq_step=step,
    )

    normalisations = miller_array.intensity_quasi_normalisations()
    normalised_intensities = miller_array.customized_copy(
        data=(miller_array.data() / normalisations.data())
    )
    reflection_table["Esq"] = flex.double(reflection_table.size(), 0.0)
    reflection_table["Esq"].set_selected(good_refl_sel, normalised_intensities.data())
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
