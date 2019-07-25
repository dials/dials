"""Test the use of BlockCalculator in ReflectionManager as used in preparation
for scan-varying refinement. This exercises the issue originally flagged in
https://github.com/dials/dials/issues/511"""

from __future__ import absolute_import, division, print_function

import pytest
from math import pi
from dials.algorithms.refinement.reflection_manager import BlockCalculator
from dxtbx.model.experiment_list import ExperimentList, Experiment
from scitbx.array_family import flex


def create_experiments(image_start=1):

    # Create models
    from libtbx.phil import parse

    overrides = """geometry.parameters.crystal.a.length.range = 10 50
  geometry.parameters.crystal.b.length.range = 10 50
  geometry.parameters.crystal.c.length.range = 10 50"""
    master_phil = parse(
        """
      include scope dials.test.algorithms.refinement.geometry_phil
      """,
        process_includes=True,
    )
    from dials.test.algorithms.refinement.setup_geometry import Extract

    models = Extract(master_phil, overrides)

    detector = models.detector
    goniometer = models.goniometer
    crystal = models.crystal
    beam = models.beam

    # Build a mock scan for a 72 degree sweep
    from dxtbx.model import ScanFactory

    sf = ScanFactory()
    scan = sf.make_scan(
        image_range=(image_start, image_start + 720 - 1),
        exposure_times=0.1,
        oscillation=(0, 0.1),
        epochs=list(range(720)),
        deg=True,
    )

    # No matter what image_start is, scan should start at 0.0 and end at 72.0 deg
    assert scan.get_oscillation_range(deg=True) == (0.0, 72.0)

    # Create an ExperimentList
    experiments = ExperimentList()
    experiments.append(
        Experiment(
            beam=beam,
            detector=detector,
            goniometer=goniometer,
            scan=scan,
            crystal=crystal,
            imageset=None,
        )
    )

    return experiments


def generate_reflections(experiments):

    from dials.algorithms.spot_prediction import IndexGenerator
    from dials.algorithms.refinement.prediction.managed_predictors import (
        ScansRayPredictor,
        ScansExperimentsPredictor,
    )
    from dials.algorithms.spot_prediction import ray_intersection
    from cctbx.sgtbx import space_group, space_group_symbols

    detector = experiments[0].detector
    crystal = experiments[0].crystal

    # All indices in a 2.0 Angstrom sphere
    resolution = 2.0
    index_generator = IndexGenerator(
        crystal.get_unit_cell(),
        space_group(space_group_symbols(1).hall()).type(),
        resolution,
    )
    indices = index_generator.to_array()

    # Predict rays within the sweep range
    scan = experiments[0].scan
    sweep_range = scan.get_oscillation_range(deg=False)
    ray_predictor = ScansRayPredictor(experiments, sweep_range)
    obs_refs = ray_predictor(indices)

    # Take only those rays that intersect the detector
    intersects = ray_intersection(detector, obs_refs)
    obs_refs = obs_refs.select(intersects)

    # Make a reflection predictor and re-predict for all these reflections. The
    # result is the same, but we gain also the flags and xyzcal.px columns
    ref_predictor = ScansExperimentsPredictor(experiments)
    obs_refs["id"] = flex.int(len(obs_refs), 0)
    obs_refs = ref_predictor(obs_refs)

    # Set 'observed' centroids from the predicted ones
    obs_refs["xyzobs.mm.value"] = obs_refs["xyzcal.mm"]

    return obs_refs


def test_per_width_and_per_image_are_equivalent():

    # Scan starting at image 1
    experiments = create_experiments(1)
    reflections = generate_reflections(experiments)

    # Check scan is consistent with the reflections
    phi_obs = reflections["xyzobs.mm.value"].parts()[2] * 180.0 / pi
    z_cal = reflections["xyzcal.px"].parts()[2]
    for phi, z in zip(phi_obs, z_cal):
        z2 = experiments[0].scan.get_array_index_from_angle(phi, deg=True)
        assert z == pytest.approx(z2)

    # Set blocks with per_width
    from copy import deepcopy

    block_calculator = BlockCalculator(experiments, deepcopy(reflections))
    im_width = experiments[0].scan.get_oscillation(deg=False)[1]
    r_pw = block_calculator.per_width(im_width, deg=False)

    # Set blocks with per_image
    block_calculator = BlockCalculator(experiments, deepcopy(reflections))
    r_pi = block_calculator.per_image()

    # Check block assignment is the same
    assert r_pw["block"].all_eq(r_pi["block"])
    for bc1, bc2 in zip(r_pw["block_centre"], r_pi["block_centre"]):
        assert bc1 == pytest.approx(bc2)

    # Scan starting at image 100
    experiments = create_experiments(100)
    reflections100 = generate_reflections(experiments)

    # Check reflections and experiments are as expected
    assert len(reflections100) == len(reflections)
    for a, b in zip(reflections, reflections100):
        assert a["xyzcal.mm"] == b["xyzcal.mm"]
    assert experiments[0].scan.get_oscillation(deg=False)[1] == im_width
    reflections = reflections100

    # Check scan is consistent with the reflections
    phi_obs = reflections["xyzobs.mm.value"].parts()[2] * 180.0 / pi
    z_cal = reflections["xyzcal.px"].parts()[2]
    for phi, z in zip(phi_obs, z_cal):
        z2 = experiments[0].scan.get_array_index_from_angle(phi, deg=True)
        assert z == pytest.approx(z2)

    # Set blocks with per_width
    block_calculator = BlockCalculator(experiments, deepcopy(reflections))
    assert experiments[0].scan.get_oscillation(deg=False)[1] == im_width
    r_pw_ = block_calculator.per_width(im_width, deg=False)

    # Block centres should have all increased by 99.0
    for a, b in zip(r_pw["block_centre"], r_pw_["block_centre"]):
        assert b == pytest.approx(a + 99.0)
    r_pw = r_pw_

    # Set blocks with per_image
    block_calculator = BlockCalculator(experiments, deepcopy(reflections))
    r_pi = block_calculator.per_image()

    # Should still give the same results as per_width
    assert r_pw["block"].all_eq(r_pi["block"])
    for bc1, bc2 in zip(r_pw["block_centre"], r_pi["block_centre"]):
        assert bc1 == pytest.approx(bc2)
