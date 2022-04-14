from __future__ import annotations

from sys import argv

import cctbx.array_family.flex
from dxtbx.model import ExperimentList

from dials.algorithms.integration.report import IntegrationReport, ProfileModelReport
from dials.algorithms.profile_model.gaussian_rs import GaussianRSProfileModeller
from dials.algorithms.profile_model.gaussian_rs import Model as GaussianRSProfileModel
from dials.algorithms.profile_model.gaussian_rs.calculator import (
    ComputeEsdBeamDivergence,
    ComputeEsdReflectingRange,
)
from dials.algorithms.shoebox import MaskCode
from dials.array_family import flex
from dials.command_line.integrate import filter_reference_pixels, process_reference
from dials.model.data import make_image
from dials.util.phil import parse
from dials_algorithms_integration_integrator_ext import ShoeboxProcessor
from dials_array_family_flex_ext import reflection_table

"""
Tested with LCY_15_15_2_*.cbf.gz using default parameters for each step of the workflow.

Kabsch 2010 refers to
Kabsch W., Integration, scaling, space-group assignment and
post-refinment, Acta Crystallographica Section D, 2010, D66, 133-144

Usage:
$ python simple_integrate.py refined.expt refined.refl
"""


if __name__ == "__main__":

    phil_scope = parse(
        """
    """
    )
    params = phil_scope.extract()

    """
    Load experiment and reflections
    """

    experiment_file = argv[1]
    reflections_file = argv[2]
    experiments = ExperimentList.from_file(experiment_file)
    reflections = reflection_table.from_msgpack_file(reflections_file)
    reflections["id"] = cctbx.array_family.flex.int(len(reflections), 0)
    reflections["imageset_id"] = cctbx.array_family.flex.int(len(reflections), 0)
    experiment = experiments[0]

    # Remove bad reflections (e.g. those not indexed)
    reflections, _ = process_reference(reflections)
    # Mask neighbouring pixels to shoeboxes
    reflections = filter_reference_pixels(reflections, experiments)

    """
    Predict reflections using experiment crystal
    """

    predicted_reflections = flex.reflection_table.from_predictions(
        experiment, padding=1.0
    )
    predicted_reflections["id"] = cctbx.array_family.flex.int(
        len(predicted_reflections), 0
    )
    predicted_reflections["imageset_id"] = cctbx.array_family.flex.int(
        len(predicted_reflections), 0
    )
    # Updates flags to set which reflections to use in generating reference profiles
    matched, reflections, unmatched = predicted_reflections.match_with_reference(
        reflections
    )

    """
    Create profile model and add it to experiment.
    This is used to predict reflection properties.
    """

    # Filter reflections to use to create the model
    min_zeta = 0.05
    used_in_ref = reflections.get_flags(reflections.flags.used_in_refinement)
    model_reflections = reflections.select(used_in_ref)
    zeta = model_reflections.compute_zeta(experiment)
    model_reflections = model_reflections.select(flex.abs(zeta) >= min_zeta)

    # sigma_D in 3.1 of Kabsch 2010
    sigma_b = ComputeEsdBeamDivergence(
        experiment.detector, model_reflections, centroid_definition="s1"
    ).sigma()

    # sigma_m in 3.1 of Kabsch 2010
    sigma_m = ComputeEsdReflectingRange(
        experiment.crystal,
        experiment.beam,
        experiment.detector,
        experiment.goniometer,
        experiment.scan,
        model_reflections,
        algorithm="extended",
    ).sigma()

    # The Gaussian model given in 2.3 of Kabsch 2010
    experiment.profile = GaussianRSProfileModel(
        params=params, n_sigma=3, sigma_b=sigma_b, sigma_m=sigma_m
    )

    """
    Compute properties for predicted reflections using profile model,
    accessed via experiment.profile_model. These reflection_table
    methods are largely just wrappers for profile_model.compute_bbox etc.

    Note: I do not think all these properties are needed for integration,
    but are all present in the current dials.integrate output.
    """

    predicted_reflections.compute_bbox(experiments)
    predicted_reflections.compute_d(experiments)
    min_zeta = 0.05
    zeta = predicted_reflections.compute_zeta(experiment)
    predicted_reflections = predicted_reflections.select(flex.abs(zeta) >= min_zeta)
    predicted_reflections.compute_partiality(experiments)

    # Shoeboxes
    predicted_reflections["shoebox"] = flex.shoebox(
        predicted_reflections["panel"],
        predicted_reflections["bbox"],
        allocate=False,
        flatten=False,
    )
    # Get actual shoebox values and the reflections for each image
    shoebox_processor = ShoeboxProcessor(
        predicted_reflections,
        len(experiment.detector),
        0,
        len(experiment.imageset),
        False,
    )

    for i in range(len(experiment.imageset)):
        image = experiment.imageset.get_corrected_data(i)
        mask = experiment.imageset.get_mask(i)
        shoebox_processor.next_data_only(make_image(image, mask))

    predicted_reflections.is_overloaded(experiments)
    predicted_reflections.compute_mask(experiments)
    predicted_reflections.contains_invalid_pixels()
    predicted_reflections.compute_background(experiments)
    predicted_reflections.compute_centroid(experiments)
    predicted_reflections.compute_summed_intensity()

    # Filter reflections with a high fraction of masked foreground
    valid_foreground_threshold = 0.75  # DIALS default
    sboxs = predicted_reflections["shoebox"]
    nvalfg = sboxs.count_mask_values(MaskCode.Valid | MaskCode.Foreground)
    nforeg = sboxs.count_mask_values(MaskCode.Foreground)
    fraction_valid = nvalfg.as_double() / nforeg.as_double()
    selection = fraction_valid < valid_foreground_threshold
    predicted_reflections.set_flags(
        selection, predicted_reflections.flags.dont_integrate
    )

    predicted_reflections["num_pixels.valid"] = sboxs.count_mask_values(MaskCode.Valid)
    predicted_reflections["num_pixels.background"] = sboxs.count_mask_values(
        MaskCode.Valid | MaskCode.Background
    )
    predicted_reflections["num_pixels.background_used"] = sboxs.count_mask_values(
        MaskCode.Valid | MaskCode.Background | MaskCode.BackgroundUsed
    )
    predicted_reflections["num_pixels.foreground"] = nvalfg

    """
    Load modeller that will calculate reference profiles and
    do the actual profile fitting integration.
    """

    # Default params when running dials.integrate with C2sum_1_*.cbf.gz
    fit_method = 1  # reciprocal space fitter (called explicitly below)
    grid_method = 2  # regular grid
    grid_size = 5  # Downsampling grid size described in 3.3 of Kabsch 2010
    num_scan_points = 72
    n_sigma = 4.5  # multiplier to expand bounding boxes
    fitting_threshold = 0.02
    reference_profile_modeller = GaussianRSProfileModeller(
        experiment.beam,
        experiment.detector,
        experiment.goniometer,
        experiment.scan,
        sigma_b,
        sigma_m,
        n_sigma,
        grid_size,
        num_scan_points,
        fitting_threshold,
        grid_method,
        fit_method,
    )

    """
    Calculate grid of reference profiles from predicted reflections
    that matched observed.
    ("Learning phase" of 3.3 in Kabsch 2010)
    """

    sel = predicted_reflections.get_flags(predicted_reflections.flags.reference_spot)
    reference_reflections = predicted_reflections.select(sel)
    sel = reference_reflections.get_flags(reference_reflections.flags.dont_integrate)
    sel = ~sel
    reference_reflections = reference_reflections.select(sel)
    reference_profile_modeller.model(reference_reflections)
    reference_profile_modeller.normalize_profiles()

    profile_model_report = ProfileModelReport(
        experiments, [reference_profile_modeller], model_reflections
    )
    print(profile_model_report.as_str(prefix=" "))

    """
    Carry out the integration by fitting to reference profiles in 1D.
    (Calculates intensity using 3.4 of Kabsch 2010)
    """

    sel = predicted_reflections.get_flags(predicted_reflections.flags.dont_integrate)
    sel = ~sel
    predicted_reflections = predicted_reflections.select(sel)
    reference_profile_modeller.fit_reciprocal_space(predicted_reflections)
    predicted_reflections.compute_corrections(experiments)

    integration_report = IntegrationReport(experiments, predicted_reflections)
    print(integration_report.as_str(prefix=" "))

    """
    Save the reflections
    """

    del predicted_reflections["shoebox"]
    sel = predicted_reflections.get_flags(
        predicted_reflections.flags.integrated, all=False
    )
    predicted_reflections = predicted_reflections.select(sel)

    predicted_reflections.as_msgpack_file("integrated.refl")
