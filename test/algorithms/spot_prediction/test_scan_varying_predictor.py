"""
Test prediction of reflections using the scan-varying reflection
predictor.
"""

from __future__ import absolute_import, division, print_function

import math

from scitbx.array_family import flex
from libtbx.phil import parse
from libtbx.test_utils import approx_equal

# Get modules to build models and minimiser using PHIL
from dials.test.algorithms.refinement import setup_geometry

# We will set up a mock scan and a mock experiment list
from dxtbx.model import ScanFactory
from dxtbx.model.experiment_list import ExperimentList, Experiment

# Reflection prediction
from dials.algorithms.spot_prediction import IndexGenerator
from dials.algorithms.spot_prediction import ray_intersection
from dials.algorithms.refinement.prediction.managed_predictors import (
    ScansRayPredictor,
    ScansExperimentsPredictor,
)
from cctbx.sgtbx import space_group, space_group_symbols


def setup_models(args):
    """setup the experimental models"""

    # Setup experimental models
    master_phil = parse(
        """
      include scope dials.test.algorithms.refinement.geometry_phil
      """,
        process_includes=True,
    )

    models = setup_geometry.Extract(master_phil, cmdline_args=args)

    detector = models.detector
    goniometer = models.goniometer
    crystal = models.crystal
    beam = models.beam

    # Build a mock scan for a 180 degree sweep
    sf = ScanFactory()
    scan = sf.make_scan(
        image_range=(1, 180),
        exposure_times=0.1,
        oscillation=(0, 1.0),
        epochs=list(range(180)),
        deg=True,
    )
    sweep_range = scan.get_oscillation_range(deg=False)
    im_width = scan.get_oscillation(deg=False)[1]
    assert sweep_range == (0.0, math.pi)
    assert approx_equal(im_width, 1.0 * math.pi / 180.0)

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


def ref_gen_static(experiments):
    """Generate some reflections using the static predictor"""

    beam = experiments[0].beam
    crystal = experiments[0].crystal
    detector = experiments[0].detector
    scan = experiments[0].scan

    # All indices to the detector max resolution
    dmin = detector.get_max_resolution(beam.get_s0())
    index_generator = IndexGenerator(
        crystal.get_unit_cell(), space_group(space_group_symbols(1).hall()).type(), dmin
    )
    indices = index_generator.to_array()

    # Predict rays within the sweep range
    sweep_range = scan.get_oscillation_range(deg=False)
    ray_predictor = ScansRayPredictor(experiments, sweep_range)
    refs = ray_predictor(indices)

    # Take only those rays that intersect the detector
    intersects = ray_intersection(detector, refs)
    refs = refs.select(intersects)

    # Make a reflection predictor and re-predict for these reflections. The
    # result is the same, but we gain also the flags and xyzcal.px columns
    ref_predictor = ScansExperimentsPredictor(experiments)
    refs["id"] = flex.int(len(refs), 0)
    refs = ref_predictor(refs)

    return refs


def ref_gen_varying(experiments):
    """Generate some reflections using the scan varying predictor"""

    beam = experiments[0].beam
    crystal = experiments[0].crystal
    detector = experiments[0].detector
    scan = experiments[0].scan

    # We need a UB matrix at the beginning of every image, and at the end of the
    # last image. These are all the same - we want to compare the scan-varying
    # predictor with the scan-static one for a flat scan.
    ar_range = scan.get_array_range()
    UBlist = [crystal.get_A() for t in range(ar_range[0], ar_range[1] + 1)]
    dmin = detector.get_max_resolution(beam.get_s0())

    from dials.algorithms.spot_prediction import ScanVaryingReflectionPredictor

    sv_predictor = ScanVaryingReflectionPredictor(experiments[0], dmin=dmin)
    refs = sv_predictor.for_ub(flex.mat3_double(UBlist))

    return refs


def sort_refs(reflections):
    """Sort reflections by Miller index and entering flag"""

    refs_sorted = sorted(reflections, key=lambda x: x["entering"])
    refs_sorted = sorted(refs_sorted, key=lambda x: x["miller_index"][2])
    refs_sorted = sorted(refs_sorted, key=lambda x: x["miller_index"][1])
    refs_sorted = sorted(refs_sorted, key=lambda x: x["miller_index"][0])

    return refs_sorted


def test():
    experiments = setup_models([])
    refs1 = ref_gen_static(experiments)
    refs2 = ref_gen_varying(experiments)

    refs1_sorted = sort_refs(refs1)
    refs2_sorted = sort_refs(refs2)

    assert len(refs1_sorted) == len(refs2_sorted)

    for (r1, r2) in zip(refs1_sorted, refs2_sorted):
        assert r1["miller_index"] == r2["miller_index"]
        dz = r1["xyzcal.px"][2] - r2["xyzcal.px"][2]
        assert abs(dz) < 0.01
