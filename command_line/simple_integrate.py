from __future__ import annotations

from sys import argv

import cctbx.array_family.flex
from dxtbx.model import ExperimentList

from dials.algorithms.profile_model.gaussian_rs import GaussianRSProfileModeller
from dials.algorithms.profile_model.gaussian_rs import Model as GaussianRSProfileModel
from dials.algorithms.profile_model.gaussian_rs.calculator import (
    ComputeEsdBeamDivergence,
    ComputeEsdReflectingRange,
)
from dials.array_family import flex
from dials.model.data import make_image
from dials.util.phil import parse
from dials_algorithms_integration_integrator_ext import ShoeboxProcessor
from dials_array_family_flex_ext import reflection_table

"""
Tested with C2sum_1_*.cbf.gz using default parameters for each step of the workflow.

Kabasch 2010 refers to
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
    experiment = experiments[0]

    """
    Predict reflections using experiment crystal
    """

    predicted_reflections = flex.reflection_table.from_predictions(experiment)
    predicted_reflections["id"] = cctbx.array_family.flex.int(
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

    reflections.compute_zeta(experiment)
    # sigma_D in 3.1 of Kabsch 2010
    sigma_b = ComputeEsdBeamDivergence(
        experiment.detector, reflections, centroid_definition="s1"
    ).sigma()
    # sigma_m in 3.1 of Kabsch 2010
    sigma_m = ComputeEsdReflectingRange(
        experiment.crystal,
        experiment.beam,
        experiment.detector,
        experiment.goniometer,
        experiment.scan,
        reflections,
    ).sigma()
    # The Gaussian model given in 2.3 of Kabsch 2010
    experiment.profile = GaussianRSProfileModel(
        params=params, n_sigma=1, sigma_b=sigma_b, sigma_m=sigma_m
    )

    """
    Compute properties for predicted reflections using profile model,
    accessed via experiments[0].profile_model. These reflection_table
    methods are just wrappers for profile_model.compute_bbox etc.
    """

    predicted_reflections.compute_bbox(experiments)
    # Shoeboxes
    predicted_reflections["shoebox"] = flex.shoebox(
        predicted_reflections["panel"],
        predicted_reflections["bbox"],
        allocate=False,
        flatten=False,
    )
    # Actual shoebox values
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

    predicted_reflections.compute_background(experiments)
    predicted_reflections.compute_centroid(experiments)
    predicted_reflections.compute_mask(experiments)
    predicted_reflections.compute_summed_intensity()
    predicted_reflections.compute_partiality(experiments)
    predicted_reflections.compute_zeta(experiment)
    predicted_reflections.compute_d_single(experiment)

    """
    Load modeller that will calculate reference profiles and
    do the actual profile fitting integration.
    """

    # Default params when running dials.integrate with C2sum_1_*.cbf.gz
    fit_method = 1  # reciprocal space fitter (called explicitly below)
    grid_method = 2  # regular grid
    grid_size = 5  # Downsampling grid size described in 3.3 of Kabsch 2010
    num_scan_points = 72
    n_sigma = 3.0  # multiplier to expand bounding boxes
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

    # Just flag all predicted reflections for integration
    predicted_reflections.unset_flags(
        flex.bool(len(predicted_reflections), True),
        predicted_reflections.flags.dont_integrate,
    )

    """
    Calculate grid of reference profiles from predicted reflections
    that matched observed.
    ("Learning phase" of 3.3 in Kabsch 2010)
    """

    reference_profile_modeller.model(predicted_reflections)
    reference_profile_modeller.finalize()

    """
    Carry out the integration by fitting to reference profiles in 1D.
    (Calculates intensity using 3.4 of Kabsch 2010)
    """

    reference_profile_modeller.fit_reciprocal_space(predicted_reflections)

    predicted_reflections.as_msgpack_file("integrated.refl")
