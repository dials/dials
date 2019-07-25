from __future__ import absolute_import, division, print_function

import sys
import pytest

from math import pi
from scitbx.array_family import flex
from dxtbx.model.experiment_list import ExperimentList, Experiment
from dials.algorithms.refinement.prediction.managed_predictors import (
    ScansRayPredictor,
    ScansExperimentsPredictor,
)
from dials.algorithms.refinement.parameterisation.scan_varying_prediction_parameters import (
    ScanVaryingPredictionParameterisation,
)
from dials.algorithms.refinement.parameterisation.scan_varying_crystal_parameters import (
    ScanVaryingCrystalOrientationParameterisation,
    ScanVaryingCrystalUnitCellParameterisation,
)
from dials.algorithms.refinement.parameterisation.scan_varying_beam_parameters import (
    ScanVaryingBeamParameterisation,
)
from dials.algorithms.refinement.parameterisation.scan_varying_detector_parameters import (
    ScanVaryingDetectorParameterisationSinglePanel,
)
from dials.algorithms.refinement.parameterisation.scan_varying_goniometer_parameters import (
    ScanVaryingGoniometerParameterisation,
)


class _Test(object):
    def create_models(self, cmdline_overrides=None):
        from dials.test.algorithms.refinement.setup_geometry import Extract
        from dxtbx.model import ScanFactory
        from libtbx.phil import parse

        if cmdline_overrides is None:
            cmdline_overrides = []
        overrides = """geometry.parameters.crystal.a.length.range = 10 50
geometry.parameters.crystal.b.length.range = 10 50
geometry.parameters.crystal.c.length.range = 10 50"""

        master_phil = parse(
            """
    include scope dials.test.algorithms.refinement.geometry_phil
    """,
            process_includes=True,
        )

        # Extract models
        models = Extract(master_phil, overrides, cmdline_args=cmdline_overrides)
        self.detector = models.detector
        self.goniometer = models.goniometer
        self.crystal = models.crystal
        self.beam = models.beam

        # Make a scan of 1-20 * 0.5 deg images
        sf = ScanFactory()
        self.scan = sf.make_scan((1, 20), 0.5, (0, 0.5), list(range(20)))

        # Generate an ExperimentList
        self.experiments = ExperimentList()
        self.experiments.append(
            Experiment(
                beam=self.beam,
                detector=self.detector,
                goniometer=self.goniometer,
                scan=self.scan,
                crystal=self.crystal,
                imageset=None,
            )
        )

        # Create a reflection predictor for the experiments
        self.ref_predictor = ScansExperimentsPredictor(self.experiments)

        # Create scan-varying parameterisations of these models, with 3 samples
        self.det_param = ScanVaryingDetectorParameterisationSinglePanel(
            self.detector, self.scan.get_array_range(), 3
        )
        self.s0_param = ScanVaryingBeamParameterisation(
            self.beam, self.scan.get_array_range(), 3, self.goniometer
        )
        self.xlo_param = ScanVaryingCrystalOrientationParameterisation(
            self.crystal, self.scan.get_array_range(), 3
        )
        self.xluc_param = ScanVaryingCrystalUnitCellParameterisation(
            self.crystal, self.scan.get_array_range(), 3
        )
        self.gon_param = ScanVaryingGoniometerParameterisation(
            self.goniometer, self.scan.get_array_range(), 3, self.beam
        )

    def generate_reflections(self):
        from cctbx.sgtbx import space_group, space_group_symbols
        from dials.algorithms.spot_prediction import IndexGenerator, ray_intersection

        sweep_range = self.scan.get_oscillation_range(deg=False)
        resolution = 2.0
        index_generator = IndexGenerator(
            self.crystal.get_unit_cell(),
            space_group(space_group_symbols(1).hall()).type(),
            resolution,
        )
        indices = index_generator.to_array()

        # Predict rays within the sweep range
        ray_predictor = ScansRayPredictor(self.experiments, sweep_range)
        obs_refs = ray_predictor(indices)

        # Take only those rays that intersect the detector
        intersects = ray_intersection(self.detector, obs_refs)
        obs_refs = obs_refs.select(intersects)

        # Re-predict using the Experiments predictor for all these reflections. The
        # result is the same, but we gain also the flags and xyzcal.px columns
        obs_refs["id"] = flex.int(len(obs_refs), 0)
        obs_refs = self.ref_predictor(obs_refs)

        # Set 'observed' centroids from the predicted ones
        obs_refs["xyzobs.mm.value"] = obs_refs["xyzcal.mm"]

        # Invent some variances for the centroid positions of the simulated data
        im_width = 0.1 * pi / 180.0
        px_size = self.detector[0].get_pixel_size()
        var_x = flex.double(len(obs_refs), (px_size[0] / 2.0) ** 2)
        var_y = flex.double(len(obs_refs), (px_size[1] / 2.0) ** 2)
        var_phi = flex.double(len(obs_refs), (im_width / 2.0) ** 2)
        obs_refs["xyzobs.mm.variance"] = flex.vec3_double(var_x, var_y, var_phi)

        # set the flex random seed to an 'uninteresting' number
        flex.set_random_seed(12407)

        # take 10 random reflections for speed
        reflections = obs_refs.select(flex.random_selection(len(obs_refs), 10))

        # use a BlockCalculator to calculate the blocks per image
        from dials.algorithms.refinement.reflection_manager import BlockCalculator

        block_calculator = BlockCalculator(self.experiments, reflections)
        reflections = block_calculator.per_image()

        return reflections


def test(cmdline_overrides=[]):
    tc = _Test()
    tc.create_models(cmdline_overrides)
    reflections = tc.generate_reflections()

    # use a ReflectionManager to exclude reflections too close to the spindle,
    # plus set the frame numbers
    from dials.algorithms.refinement.reflection_manager import ReflectionManager

    refman = ReflectionManager(reflections, tc.experiments, outlier_detector=None)
    refman.finalise()

    # create prediction parameterisation of the requested type
    pred_param = ScanVaryingPredictionParameterisation(
        tc.experiments,
        [tc.det_param],
        [tc.s0_param],
        [tc.xlo_param],
        [tc.xluc_param],
        [tc.gon_param],
    )

    # keep only those reflections that pass inclusion criteria and have predictions
    reflections = refman.get_matches()

    # get analytical gradients
    pred_param.compose(reflections)
    an_grads = pred_param.get_gradients(reflections)

    # get finite difference gradients
    p_vals = pred_param.get_param_vals()
    p_names = pred_param.get_param_names()
    deltas = [1.0e-7] * len(p_vals)

    for i, delta in enumerate(deltas):
        val = p_vals[i]

        p_vals[i] -= delta / 2.0
        pred_param.set_param_vals(p_vals)
        pred_param.compose(reflections)

        tc.ref_predictor(reflections)

        rev_state = reflections["xyzcal.mm"].deep_copy()

        p_vals[i] += delta
        pred_param.set_param_vals(p_vals)
        pred_param.compose(reflections)

        tc.ref_predictor(reflections)

        fwd_state = reflections["xyzcal.mm"].deep_copy()
        p_vals[i] = val

        fd = fwd_state - rev_state
        x_grads, y_grads, phi_grads = fd.parts()
        x_grads /= delta
        y_grads /= delta
        phi_grads /= delta

        try:
            for (a, b) in zip(x_grads, an_grads[i]["dX_dp"]):
                assert a == pytest.approx(b, abs=1e-5)
            for (a, b) in zip(y_grads, an_grads[i]["dY_dp"]):
                assert a == pytest.approx(b, abs=1e-5)
            for (a, b) in zip(phi_grads, an_grads[i]["dphi_dp"]):
                assert a == pytest.approx(b, abs=1e-5)
        except AssertionError:
            print("Failure for {}".format(p_names[i]))
            raise

    # return to the initial state
    pred_param.set_param_vals(p_vals)
    pred_param.compose(reflections)


if __name__ == "__main__":
    cmdline_overrides = sys.argv[1:]
    test(cmdline_overrides)
