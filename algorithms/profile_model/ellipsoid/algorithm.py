#!/usr/bin/env python
#
# dials.ellipsoid.py
#
#  Copyright (C) 2018 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import annotations

import copy
import logging
from math import pi

import numpy as np
from numpy.linalg import inv, norm

from dxtbx import flumpy
from libtbx.phil import parse
from scitbx import matrix

from dials.algorithms.profile_model.ellipsoid import chisq_pdf
from dials.algorithms.profile_model.ellipsoid.model import ProfileModelFactory
from dials.algorithms.profile_model.ellipsoid.parameterisation import ModelState
from dials.algorithms.profile_model.ellipsoid.refiner import Refiner as ProfileRefiner
from dials.algorithms.profile_model.ellipsoid.refiner import RefinerData
from dials.algorithms.profile_model.gaussian_rs import BBoxCalculator, MaskCalculator
from dials.algorithms.profile_model.gaussian_rs.calculator import (
    ComputeEsdBeamDivergence,
)
from dials.algorithms.spot_prediction import IndexGenerator
from dials.array_family import flex

logger = logging.getLogger("dials")

ellipsoid_algorithm_phil_scope = parse(
    """
  integration {

    use_crude_shoebox_mask = False
      .type = bool

    shoebox {

      probability = 0.9973
        .type = float

    }

    corrections {
      lp = True
        .type = bool
      dqe = True
        .type = bool
      partiality = True
        .type = bool
    }

  }

  debug {
    output {

      strong_spots = False
        .type = bool

      shoeboxes = False
        .type = bool

      profile_model = True
        .type = bool

      history = True
        .type = bool

      plots = False
        .type = bool

      print_shoeboxes = False
        .type = bool
    }
  }


"""
)

# Parameters
phil_scope = parse(
    """
include scope dials.algorithms.profile_model.ellipsoid.model.phil_scope
""",
    process_includes=True,
).adopt_scope(ellipsoid_algorithm_phil_scope)


def _compute_beam_vector(experiment, reflection_table):
    """Compute the obseved beam vector"""
    s1_obs = flex.vec3_double(len(reflection_table))
    for i in range(len(s1_obs)):
        x, y, _ = reflection_table["xyzobs.px.value"][i]
        s1_obs[i] = experiment.detector[0].get_pixel_lab_coord((x, y))
    return s1_obs


def _compute_bbox(experiment, reflection_table, sigma_d, s1="s1_obs"):
    """Compute the bounding box"""

    # Initialise the bounding box calculator
    compute_bbox = BBoxCalculator(
        experiment.crystal,
        experiment.beam,
        experiment.detector,
        experiment.goniometer,
        experiment.scan,
        sigma_d * 6,
        0,
    )

    # Compute the bounding box
    bbox = compute_bbox(
        reflection_table[s1],
        reflection_table["xyzcal.px"].parts()[2],
        reflection_table["panel"],
    )

    return bbox


def _compute_mask(
    experiment, reflection_table, sigma_d, s1="s1_obs", strong_shoeboxes=None
):
    """Compute the spot mask"""

    # Initialise the mask calculator
    mask_foreground = MaskCalculator(
        experiment.crystal,
        experiment.beam,
        experiment.detector,
        experiment.goniometer,
        experiment.scan,
        sigma_d * 3,
        0,
    )

    # Compute the reflection mask
    mask_foreground(
        reflection_table["shoebox"],
        reflection_table[s1],
        reflection_table["xyzcal.px"].parts()[2],
        reflection_table["panel"],
    )

    # if strong_shoeboxes given, apply the mask from this
    if strong_shoeboxes:
        # Apply strong spot mask
        assert len(reflection_table) == len(strong_shoeboxes)
        new_shoeboxes = reflection_table["shoebox"]
        for s in range(len(new_shoeboxes)):
            bbox_old = strong_shoeboxes[s].bbox
            mask_old = strong_shoeboxes[s].mask
            bbox_new = new_shoeboxes[s].bbox
            mask_new = new_shoeboxes[s].mask
            for j in range(mask_old.all()[1]):
                for i in range(mask_old.all()[2]):
                    ii = bbox_old[0] + i - bbox_new[0]
                    jj = bbox_old[2] + j - bbox_new[2]
                    if mask_old[0, j, i] == 5:
                        if (
                            ii >= 0
                            and jj >= 0
                            and jj < mask_new.all()[1]
                            and ii < mask_new.all()[2]
                        ):
                            assert mask_new[0, jj, ii] & (1 << 0)
                            mask_new[0, jj, ii] |= (1 << 2) | (1 << 3)


def initial_integrator(experiments, reflection_table):
    """Performs an initial integration of strong spots"""

    # some functions require an experimentlist, others just the experiment
    experiment = experiments[0]
    sel = reflection_table.get_flags(reflection_table.flags.strong)
    strong_refls = reflection_table.select(sel)
    strong_shoeboxes = copy.deepcopy(
        strong_refls["shoebox"]
    )  # Save the strong shoeboxes

    # Compute and initial spot size estimate and beam vector
    sigma_d = ComputeEsdBeamDivergence(experiment.detector, strong_refls).sigma()
    sigma_degrees = sigma_d * 180 / pi
    logger.info(
        f"Initial sigma d estimate for {len(strong_refls)} reflections\n"
        + f"Sigma D: {sigma_degrees:.5f} degrees\n",
    )
    strong_refls["s1_obs"] = _compute_beam_vector(experiment, strong_refls)
    strong_refls["bbox"] = _compute_bbox(experiment, strong_refls, sigma_d, "s1_obs")

    # allocate and extract shoebox
    strong_refls["shoebox"] = flex.shoebox(
        strong_refls["panel"], strong_refls["bbox"], allocate=True
    )
    strong_refls.extract_shoeboxes(experiment.imageset)

    _compute_mask(experiment, strong_refls, sigma_d, "s1_obs", strong_shoeboxes)

    logger.info(
        f"Computing background, intensity and centroids for {len(strong_refls)} reflections"
    )
    strong_refls.compute_background(experiments)
    strong_refls.compute_summed_intensity()

    n_sum = strong_refls.get_flags(strong_refls.flags.integrated_sum).count(True)
    logger.info(f"{n_sum} reflections integrated")

    strong_refls.compute_centroid(experiments)
    strong_refls.compute_d(experiments)

    return strong_refls, sigma_d


def final_integrator(
    experiments,
    reflection_table,
    sigma_d,
    use_crude_shoebox_mask=False,
    shoebox_probability=0.9973,
):
    """Performs an initial integration of all predicted spots"""

    experiment = experiments[0]
    logger.info("\n" + "=" * 80 + "\nIntegrating reflections")

    # first compute the bbox
    if use_crude_shoebox_mask:
        reflection_table["bbox"] = _compute_bbox(
            experiment, reflection_table, sigma_d, "s1"
        )
    else:
        # compute bbox from model
        profile = experiment.crystal.mosaicity
        profile.parameterisation.compute_bbox(
            experiments, reflection_table, shoebox_probability
        )

    # Select reflections within detector
    x0, x1, y0, y1, _, _ = reflection_table["bbox"].parts()
    xsize, ysize = experiment.detector[0].get_image_size()
    selection = (x1 > 0) & (y1 > 0) & (x0 < xsize) & (y0 < ysize)
    reflection_table = reflection_table.select(selection)

    reflection_table["shoebox"] = flex.shoebox(
        reflection_table["panel"], reflection_table["bbox"], allocate=True
    )
    reflection_table.extract_shoeboxes(experiment.imageset)

    if use_crude_shoebox_mask:
        _compute_mask(experiment, reflection_table, sigma_d, "s1")
    else:
        # compute the mask from the model.
        profile = experiment.crystal.mosaicity
        profile.parameterisation.compute_mask(
            experiments, reflection_table, shoebox_probability
        )

    logger.info(
        f"Computing background, intensity, corrections for {len(reflection_table)} reflections"
    )
    reflection_table.compute_background(experiments)
    reflection_table.compute_summed_intensity()
    reflection_table.compute_corrections(experiments)
    reflection_table.compute_centroid(experiments)

    profile = experiment.crystal.mosaicity
    profile.parameterisation.compute_partiality(experiments, reflection_table)

    return reflection_table


def refine_profile(experiment, profile, refiner_data, wavelength_spread_model="delta"):
    """Do the profile refinement"""
    logger.info("\n" + "=" * 80 + "\nRefining profile parmameters")

    # Create the parameterisation
    state = ModelState(
        experiment,
        profile.parameterisation.parameterisation(),
        fix_orientation=True,
        fix_unit_cell=True,
        fix_wavelength_spread=wavelength_spread_model == "delta",
    )

    # Create the refiner and refine
    refiner = ProfileRefiner(state, refiner_data)
    refiner.refine()

    # Set the profile parameters
    profile.parameterisation.update_model(state)
    # Set the mosaicity
    experiment.crystal.mosaicity = profile

    return refiner


def refine_crystal(
    experiment,
    profile,
    refiner_data,
    fix_unit_cell=False,
    fix_orientation=False,
    wavelength_spread_model="delta",
):
    """Do the crystal refinement"""
    if (fix_unit_cell is True) and (fix_orientation is True):
        return

    logger.info("\n" + "=" * 80 + "\nRefining crystal parmameters")

    # Create the parameterisation
    state = ModelState(
        experiment,
        profile.parameterisation.parameterisation(),
        fix_mosaic_spread=True,
        fix_unit_cell=fix_unit_cell,
        fix_orientation=fix_orientation,
        fix_wavelength_spread=wavelength_spread_model == "delta",
    )

    # Create the refiner and refine
    refiner = ProfileRefiner(state, refiner_data)
    refiner.refine()

    return refiner


def predict_after_ellipsoid_refinement(experiment, reflection_table):
    """
    Predict the position of the spots

    """
    # Get some stuff from experiment
    A = np.array(experiment.crystal.get_A(), dtype=np.float64).reshape((3, 3))
    s0 = np.array([experiment.beam.get_s0()], dtype=np.float64).reshape(3, 1)
    s0_length = norm(s0)

    # Compute the vector to the reciprocal lattice point
    # since this is not on the ewald sphere, lets call it s2
    h = reflection_table["miller_index"]
    s1 = flex.vec3_double(len(h))
    s2 = flex.vec3_double(len(h))
    for i in range(len(reflection_table)):
        r = np.matmul(A, np.array([h[i]], dtype=np.float64).reshape(3, 1))
        s2_i = r + s0
        s2[i] = matrix.col(flumpy.from_numpy(s2_i))
        s1[i] = matrix.col(flumpy.from_numpy(s2_i * s0_length / norm(s2_i)))
    reflection_table["s1"] = s1
    reflection_table["s2"] = s2
    reflection_table["entering"] = flex.bool(len(h), False)

    # Compute the ray intersections
    xyzpx = flex.vec3_double()
    xyzmm = flex.vec3_double()
    for i in range(len(s2)):
        ss = s1[i]
        mm = experiment.detector[0].get_ray_intersection(ss)
        px = experiment.detector[0].millimeter_to_pixel(mm)
        xyzpx.append(px + (0,))
        xyzmm.append(mm + (0,))
    reflection_table["xyzcal.mm"] = xyzmm
    reflection_table["xyzcal.px"] = xyzpx
    return reflection_table


def compute_prediction_probability(experiment, reflection_table):

    # Get stuff from experiment
    s0 = np.array([experiment.beam.get_s0()], dtype=np.float64).reshape(3, 1)
    s0_length = norm(s0)
    profile = experiment.crystal.mosaicity

    # Loop through reflections
    min_p = None
    for s2 in reflection_table["s2"]:
        s2 = np.array(list(s2), dtype=np.float64).reshape(3, 1)
        s3 = s2 * s0_length / norm(s2)
        r = s2 - s0
        epsilon = s3 - s2
        sigma = profile.parameterisation.sigma_for_reflection(s0, r)
        sigma_inv = inv(sigma)
        d = np.matmul(np.matmul(epsilon.T, sigma_inv), epsilon)[0, 0]
        p = chisq_pdf(3, float(d))
        if min_p is None or p < min_p:
            min_p = p

    # Print some stuff
    logger.info(
        "Quantile required to predicted all observed reflections = %.5f" % (1 - min_p)
    )


def run_ellipsoid_refinement(
    experiments,
    reflection_table,
    sigma_d,
    profile_model="angular2",
    wavelength_model="delta",
    fix_unit_cell=False,
    fix_orientation=False,
    capture_progress=False,
    n_cycles=3,
):
    """Runs ellipsoid refinement on strong spots.

    Creates the necessary data needed, then runs cycles of profile and crystal
    refinement,"""

    output_data = {
        "refiner_output": {
            "history": [],
            "correlation": None,
            "labels": None,
        },
    }

    # Set the M params
    if not hasattr(experiments[0].crystal, "mosaicity"):
        profile = ProfileModelFactory.from_sigma_d(profile_model, sigma_d)
    else:
        profile = experiments[0].crystal.mosaicity

    # Construct the profile refiner data
    refiner_data = RefinerData.from_reflections(experiments[0], reflection_table)

    # Do the refinement
    for _ in range(n_cycles):
        # refine the profile
        refiner = refine_profile(
            experiments[0],
            profile,
            refiner_data,
            wavelength_spread_model=wavelength_model,
        )
        if capture_progress:
            # Save some data for plotting later.
            output_data["refiner_output"]["history"].append(refiner.history)
            output_data["refiner_output"]["correlation"] = refiner.correlation()
            output_data["refiner_output"]["labels"] = refiner.labels()

        # refine the crystal
        _ = refine_crystal(
            experiments[0],
            profile,
            refiner_data,
            fix_unit_cell=fix_unit_cell,
            fix_orientation=fix_orientation,
            wavelength_spread_model=wavelength_model,
        )

    experiments[0].profile = profile

    # Post process the reflections
    # Update predictions
    reflection_table = predict_after_ellipsoid_refinement(
        experiments[0], reflection_table
    )
    # Compute prob
    compute_prediction_probability(experiments[0], reflection_table)

    return experiments, reflection_table, output_data


def predict(experiments, d_min=None, prediction_probability=0.9973):
    """Predict the reflections"""
    logger.info("\n" + "=" * 80 + "\nPredicting reflections")

    # Set a resolution range
    if d_min is None:
        s0 = experiments[0].beam.get_s0()
        d_min = experiments[0].detector.get_max_resolution(s0)

    # Create the index generator
    index_generator = IndexGenerator(
        experiments[0].crystal.get_unit_cell(),
        experiments[0].crystal.get_space_group().type(),
        d_min,
    )

    # Get an array of miller indices
    miller_indices_to_test = index_generator.to_array()
    logger.info("Generated %d miller indices" % len(miller_indices_to_test))

    # Get the covariance matrix
    profile = experiments[0].crystal.mosaicity

    reflection_table = profile.parameterisation.predict_reflections(
        experiments, miller_indices_to_test, prediction_probability
    )

    # Do the prediction
    reflection_table.compute_d(experiments)
    logger.info("Predicted %d reflections" % len(reflection_table))

    return reflection_table
