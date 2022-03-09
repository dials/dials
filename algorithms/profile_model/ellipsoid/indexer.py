from __future__ import annotations

import logging
from math import floor, sqrt

import numpy as np
from numpy.linalg import inv, norm

from cctbx.array_family import flex
from dxtbx import flumpy
from scitbx import matrix

from dials.algorithms.profile_model.ellipsoid import chisq_quantile
from dials.algorithms.statistics.fast_mcd import FastMCD, maha_dist_sq

logger = logging.getLogger("dials")


def _index(reflection_table, experiment, fail_on_bad_index=False):
    """Index the strong spots"""

    # Get some stuff from experiment
    A = np.array(experiment.crystal.get_A(), dtype=np.float64).reshape(3, 3)
    s0 = np.array([experiment.beam.get_s0()], dtype=np.float64).reshape(3, 1)
    s0_length = norm(s0)
    detector = experiment.detector

    # Create array if necessary
    if "miller_index" not in reflection_table:
        reflection_table["miller_index"] = flex.miller_index(len(reflection_table))

    # Index all the reflections
    miller_index = reflection_table["miller_index"]
    selection = flex.size_t()
    num_reindexed = 0
    for i, xyz in enumerate(reflection_table["xyzobs.px.value"]):
        # Get the observed pixel coordinate
        x, y, _ = xyz

        # Get the lab coord
        s1 = np.array(
            detector[0].get_pixel_lab_coord((x, y)), dtype=np.float64
        ).reshape(3, 1)
        s1_norm = norm(s1)
        s1 *= s0_length / s1_norm

        # Get the reciprocal lattice vector
        r = s1 - s0
        # Compute the fractional miller index
        hf = np.matmul(inv(A), r)
        # Compute the integer miller index
        h = np.array([int(floor(j + 0.5)) for j in hf[:, 0]], dtype=int).reshape(3, 1)

        # Print warning if reindexing
        if tuple(h) != miller_index[i]:
            logger.warn(
                "Reindexing (% 3d, % 3d, % 3d) -> (% 3d, % 3d, % 3d)"
                % (miller_index[i] + tuple(h))
            )
            num_reindexed += 1
            miller_index[i] = matrix.col(flumpy.from_numpy(h))
            if fail_on_bad_index:
                raise RuntimeError("Bad index")

        # If its not indexed as 0, 0, 0 then append
        if h.any() and norm(h - hf) < 0.3:
            selection.append(i)

    # Print some info
    logger.info(
        "Reindexed %d/%d input reflections" % (num_reindexed, len(reflection_table))
    )
    logger.info(
        "Selected %d/%d input reflections" % (len(selection), len(reflection_table))
    )

    # Select all the indexed reflections
    reflection_table.set_flags(selection, reflection_table.flags.indexed)
    reflection_table = reflection_table.select(selection)
    return reflection_table


def _predict(reflection_table, experiment):
    """
    Predict the position of the spots

    """

    # Get some stuff from experiment
    A = np.array(experiment.crystal.get_A(), dtype=np.float64).reshape((3, 3))
    s0 = np.array([experiment.beam.get_s0()], dtype=np.float64).reshape(3, 1)
    s0_length = norm(s0)

    # Compute the vector to the reciprocal lattice point
    # since this is not on the ewald sphere, lets call it s2
    s1 = flex.vec3_double(reflection_table.size())
    s2 = flex.vec3_double(reflection_table.size())
    for i, h in enumerate(reflection_table["miller_index"]):
        r = np.matmul(A, np.array([h], dtype=np.float64).reshape(3, 1))
        s2_i = r + s0
        s2[i] = matrix.col(flumpy.from_numpy(s2_i))
        s1[i] = matrix.col(flumpy.from_numpy(s2_i * s0_length / norm(s2_i)))
    reflection_table["s1"] = s1
    reflection_table["s2"] = s2
    reflection_table["entering"] = flex.bool(reflection_table.size(), False)

    # Compute the ray intersections
    xyzpx = flex.vec3_double()
    xyzmm = flex.vec3_double()
    for ss in s1:
        mm = experiment.detector[0].get_ray_intersection(ss)
        px = experiment.detector[0].millimeter_to_pixel(mm)
        xyzpx.append(px + (0,))
        xyzmm.append(mm + (0,))
    reflection_table["xyzcal.mm"] = xyzmm
    reflection_table["xyzcal.px"] = xyzpx
    logger.info("Do prediction for %d reflections" % len(reflection_table))
    return reflection_table


def _filter_reflections_based_on_centroid_distance(
    reflection_table,
    experiment,
    outlier_probability=0.975,
    max_separation=2,
):
    """
    Filter reflections too far from predicted position

    """

    # Compute the x and y residuals
    Xobs, Yobs, _ = reflection_table["xyzobs.px.value"].parts()
    Xcal, Ycal, _ = reflection_table["xyzcal.px"].parts()
    Xres = Xobs - Xcal
    Yres = Yobs - Ycal

    # Compute the epsilon residual
    s0_length = 1.0 / experiment.beam.get_wavelength()
    s1x, s1y, s1z = reflection_table["s2"].parts()
    s1_length = flex.sqrt(s1x**2 + s1y**2 + s1z**2)
    Eres = s1_length - s0_length

    # Initialise the fast_mcd outlier algorithm
    # fast_mcd = FastMCD((Xres, Yres, Eres))
    fast_mcd = FastMCD((Xres, Yres))

    # get location and MCD scatter estimate
    T, S = fast_mcd.get_corrected_T_and_S()

    # get squared Mahalanobis distances
    # d2s = maha_dist_sq((Xres, Yres, Eres), T, S)
    d2s = maha_dist_sq((Xres, Yres), T, S)

    # Compute the cutoff
    mahasq_cutoff = chisq_quantile(2, outlier_probability)

    # compare to the threshold and select reflections
    selection1 = d2s < mahasq_cutoff
    selection2 = flex.sqrt(Xres**2 + Yres**2) < max_separation
    selection = selection1 & selection2
    reflection_table = reflection_table.select(selection)
    n_refl = reflection_table.size()

    # Print some stuff
    logger.info("-" * 80)
    logger.info("Centroid outlier rejection")
    logger.info(f" Using MCD algorithm with probability = {outlier_probability}")
    logger.info(" Max X residual: %f" % flex.max(flex.abs(Xres)))
    logger.info(" Max Y residual: %f" % flex.max(flex.abs(Yres)))
    logger.info(" Max E residual: %f" % flex.max(flex.abs(Eres)))
    logger.info(" Mean X RMSD: %f" % (sqrt(flex.sum(Xres**2) / len(Xres))))
    logger.info(" Mean Y RMSD: %f" % (sqrt(flex.sum(Yres**2) / len(Yres))))
    logger.info(" Mean E RMSD: %f" % (sqrt(flex.sum(Eres**2) / len(Eres))))
    logger.info(" MCD location estimate: %.4f, %.4f" % tuple(T))
    logger.info(
        """ MCD scatter estimate:
    %.7f, %.7f,
    %.7f, %.7f"""
        % tuple(S)
    )
    logger.info(" Number of outliers: %d" % selection1.count(False))
    logger.info(
        " Number of reflections with residual > %0.2f pixels: %d"
        % (max_separation, selection2.count(False))
    )
    logger.info(f"Number of reflections selection for refinement: {n_refl}")
    logger.info("-" * 80)

    return reflection_table


def reindex(
    reflection_table,
    experiment,
    outlier_probability=0.975,
    max_separation=2,
    fail_on_bad_index=False,
):
    """Reindex strong spots and perform filtering"""
    reflection_table = _index(reflection_table, experiment, fail_on_bad_index)
    reflection_table = _predict(reflection_table, experiment)
    reflection_table = _filter_reflections_based_on_centroid_distance(
        reflection_table,
        experiment,
        outlier_probability=outlier_probability,
        max_separation=max_separation,
    )
    return reflection_table
