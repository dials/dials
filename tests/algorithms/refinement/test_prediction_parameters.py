from __future__ import annotations

import math
from os.path import join

import pytest

from cctbx.sgtbx import space_group, space_group_symbols
from dxtbx.format.FormatISISSXD import FormatISISSXD
from dxtbx.model import CrystalFactory
from dxtbx.model.experiment_list import Experiment, ExperimentList
from libtbx.phil import parse
from scitbx.array_family import flex

from dials.algorithms.refinement.parameterisation.beam_parameters import (
    BeamParameterisation,
)
from dials.algorithms.refinement.parameterisation.crystal_parameters import (
    CrystalOrientationParameterisation,
    CrystalUnitCellParameterisation,
)
from dials.algorithms.refinement.parameterisation.detector_parameters import (
    DetectorParameterisationHierarchical,
    DetectorParameterisationSinglePanel,
)
from dials.algorithms.refinement.parameterisation.goniometer_parameters import (
    GoniometerParameterisation,
)
from dials.algorithms.refinement.parameterisation.prediction_parameters import (
    LauePredictionParameterisation,
    XYPhiPredictionParameterisation,
)
from dials.algorithms.refinement.prediction.managed_predictors import (
    LaueExperimentsPredictor,
    ScansExperimentsPredictor,
    ScansRayPredictor,
)
from dials.algorithms.refinement.reflection_manager import LaueReflectionManager
from dials.algorithms.spot_prediction import (
    IndexGenerator,
    LaueReflectionPredictor,
    ray_intersection,
)

from . import geometry_phil
from .setup_geometry import Extract


def test():
    overrides = """geometry.parameters.crystal.a.length.range = 10 50
  geometry.parameters.crystal.b.length.range = 10 50
  geometry.parameters.crystal.c.length.range = 10 50"""

    master_phil = parse(geometry_phil)

    models = Extract(master_phil, overrides)

    mydetector = models.detector
    mygonio = models.goniometer
    mycrystal = models.crystal
    mybeam = models.beam

    # Build a mock scan for a 72 degree sequence
    sequence_range = (0.0, math.pi / 5.0)
    from dxtbx.model import ScanFactory

    sf = ScanFactory()
    myscan = sf.make_scan(
        image_range=(1, 720),
        exposure_times=0.1,
        oscillation=(0, 0.1),
        epochs=list(range(720)),
        deg=True,
    )

    #### Create parameterisations of these models
    det_param = DetectorParameterisationSinglePanel(mydetector)
    s0_param = BeamParameterisation(mybeam, mygonio)
    xlo_param = CrystalOrientationParameterisation(mycrystal)
    xluc_param = CrystalUnitCellParameterisation(mycrystal)
    gon_param = GoniometerParameterisation(mygonio, mybeam)

    # Create an ExperimentList
    experiments = ExperimentList()
    experiments.append(
        Experiment(
            beam=mybeam,
            detector=mydetector,
            goniometer=mygonio,
            scan=myscan,
            crystal=mycrystal,
            imageset=None,
        )
    )

    #### Unit tests

    # Build a prediction parameterisation
    pred_param = XYPhiPredictionParameterisation(
        experiments,
        detector_parameterisations=[det_param],
        beam_parameterisations=[s0_param],
        xl_orientation_parameterisations=[xlo_param],
        xl_unit_cell_parameterisations=[xluc_param],
        goniometer_parameterisations=[gon_param],
    )

    # Generate reflections
    resolution = 2.0
    index_generator = IndexGenerator(
        mycrystal.get_unit_cell(),
        space_group(space_group_symbols(1).hall()).type(),
        resolution,
    )
    indices = index_generator.to_array()

    # Predict rays within the sequence range
    ray_predictor = ScansRayPredictor(experiments, sequence_range)
    obs_refs = ray_predictor(indices)

    # Take only those rays that intersect the detector
    intersects = ray_intersection(mydetector, obs_refs)
    obs_refs = obs_refs.select(intersects)

    # Make a reflection predictor and re-predict for all these reflections. The
    # result is the same, but we gain also the flags and xyzcal.px columns
    ref_predictor = ScansExperimentsPredictor(experiments)
    obs_refs["id"] = flex.int(len(obs_refs), 0)
    obs_refs = ref_predictor(obs_refs)

    # Set 'observed' centroids from the predicted ones
    obs_refs["xyzobs.mm.value"] = obs_refs["xyzcal.mm"]

    # Invent some variances for the centroid positions of the simulated data
    im_width = 0.1 * math.pi / 180.0
    px_size = mydetector[0].get_pixel_size()
    var_x = flex.double(len(obs_refs), (px_size[0] / 2.0) ** 2)
    var_y = flex.double(len(obs_refs), (px_size[1] / 2.0) ** 2)
    var_phi = flex.double(len(obs_refs), (im_width / 2.0) ** 2)
    obs_refs["xyzobs.mm.variance"] = flex.vec3_double(var_x, var_y, var_phi)

    # use a ReflectionManager to exclude reflections too close to the spindle
    from dials.algorithms.refinement.reflection_manager import ReflectionManager

    refman = ReflectionManager(obs_refs, experiments, outlier_detector=None)
    refman.finalise()

    # Redefine the reflection predictor to use the type expected by the Target class
    ref_predictor = ScansExperimentsPredictor(experiments)

    # keep only those reflections that pass inclusion criteria and have predictions
    reflections = refman.get_matches()

    # get analytical gradients
    an_grads = pred_param.get_gradients(reflections)

    # get finite difference gradients
    p_vals = pred_param.get_param_vals()
    deltas = [1.0e-7] * len(p_vals)

    for i, delta in enumerate(deltas):
        val = p_vals[i]

        p_vals[i] -= delta / 2.0
        pred_param.set_param_vals(p_vals)

        ref_predictor(reflections)

        rev_state = reflections["xyzcal.mm"].deep_copy()

        p_vals[i] += delta
        pred_param.set_param_vals(p_vals)

        ref_predictor(reflections)

        fwd_state = reflections["xyzcal.mm"].deep_copy()
        p_vals[i] = val

        fd = fwd_state - rev_state
        x_grads, y_grads, phi_grads = fd.parts()
        x_grads /= delta
        y_grads /= delta
        phi_grads /= delta

        # compare with analytical calculation
        assert x_grads == pytest.approx(an_grads[i]["dX_dp"], abs=5.0e-6)
        assert y_grads == pytest.approx(an_grads[i]["dY_dp"], abs=5.5e-6)
        assert phi_grads == pytest.approx(an_grads[i]["dphi_dp"], abs=5.0e-6)

    # return to the initial state
    pred_param.set_param_vals(p_vals)


def test_laue_prediction_parameters(dials_data):
    fmt = FormatISISSXD(join(dials_data("isis_sxd_example_data"), "sxd_nacl_run.nxs"))
    beam = fmt.get_beam()
    detector = fmt.get_detector()
    goniometer = fmt.get_goniometer()
    scan = fmt.get_scan()
    crystal = CrystalFactory.from_dict(
        {
            "__id__": "crystal",
            "real_space_a": (
                0.5681647125795644,
                -2.9735716012061135,
                -2.707784412005687,
            ),
            "real_space_b": (
                -2.4994848902125884,
                -2.3900344014694066,
                2.091613643314567,
            ),
            "real_space_c": (
                -1.2771711635863638,
                3.676428861690809,
                -1.226011051463438,
            ),
            "space_group_hall_symbol": " P 1",
            "B_covariance": (
                2.618491627225783e-13,
                -2.4190170785778272e-30,
                2.7961382012436816e-30,
                1.4283218313839273e-13,
                8.110824693143866e-15,
                2.7961382012436816e-30,
                -1.922218398881239e-13,
                -1.1641948761717081e-14,
                2.2832201114561855e-14,
                -2.419017078577827e-30,
                1.3543505986455804e-44,
                -8.081590630292518e-46,
                -4.202632560757537e-29,
                -5.437640708903305e-29,
                -8.081590630292518e-46,
                3.330706229067803e-30,
                5.621471188408899e-29,
                -6.599119546892406e-30,
                2.7961382012436816e-30,
                -8.08159063029252e-46,
                9.550033948814972e-46,
                5.487666450546843e-30,
                2.7096475027184553e-30,
                9.550033948814972e-46,
                -3.935814660390771e-30,
                -3.889472044173952e-30,
                7.798194512461942e-30,
                1.428321831383927e-13,
                -4.2026325607575364e-29,
                5.487666450546843e-30,
                7.789867544667339e-13,
                1.4101250207277487e-13,
                5.487666450546843e-30,
                -2.0005409484272627e-13,
                -2.021584892435437e-13,
                4.481019714719027e-14,
                8.110824693143867e-15,
                -5.437640708903304e-29,
                2.7096475027184553e-30,
                1.4101250207277487e-13,
                2.5553690436147e-13,
                2.7096475027184553e-30,
                -1.1167612085554417e-14,
                -1.8848015530742402e-13,
                2.2125950964841596e-14,
                2.7961382012436816e-30,
                -8.08159063029252e-46,
                9.550033948814972e-46,
                5.487666450546843e-30,
                2.7096475027184553e-30,
                9.550033948814972e-46,
                -3.935814660390771e-30,
                -3.889472044173952e-30,
                7.798194512461942e-30,
                -1.922218398881239e-13,
                3.330706229067804e-30,
                -3.93581466039077e-30,
                -2.000540948427263e-13,
                -1.1167612085554417e-14,
                -3.93581466039077e-30,
                2.7092227778026175e-13,
                1.6029668235488112e-14,
                -3.2138365634328507e-14,
                -1.1641948761717081e-14,
                5.621471188408898e-29,
                -3.889472044173952e-30,
                -2.021584892435437e-13,
                -1.88480155307424e-13,
                -3.889472044173952e-30,
                1.6029668235488112e-14,
                2.7054780216756276e-13,
                -3.175994945548343e-14,
                2.2832201114561858e-14,
                -6.599119546892407e-30,
                7.79819451246194e-30,
                4.4810197147190265e-14,
                2.2125950964841592e-14,
                7.79819451246194e-30,
                -3.2138365634328507e-14,
                -3.175994945548343e-14,
                6.36770905528953e-14,
            ),
        }
    )

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

    det_param = DetectorParameterisationHierarchical(detector)
    xlo_param = CrystalOrientationParameterisation(crystal)
    xluc_param = CrystalUnitCellParameterisation(crystal)

    pred_param = LauePredictionParameterisation(
        experiments,
        detector_parameterisations=[det_param],
        beam_parameterisations=[],
        xl_orientation_parameterisations=[xlo_param],
        xl_unit_cell_parameterisations=[xluc_param],
    )

    reflection_predictor = LaueReflectionPredictor(experiments[0], 1.0)
    obs_refs = reflection_predictor.all_reflections_for_asu(0.0)

    # Set 'observed' centroids from the predicted ones
    obs_refs["xyzobs.mm.value"] = obs_refs["xyzcal.mm"]
    obs_refs["s0"] = obs_refs["s0_cal"]
    obs_refs["wavelength"] = obs_refs["wavelength_cal"]
    obs_refs["id"] = flex.int(len(obs_refs), 0)

    # Invent some variances for the centroid positions of the simulated data
    px_size = detector[0].get_pixel_size()
    var_x = flex.double(len(obs_refs), (px_size[0] / 2.0) ** 2)
    var_y = flex.double(len(obs_refs), (px_size[1] / 2.0) ** 2)
    var_z = flex.double(len(obs_refs), 0.0)
    obs_refs["xyzobs.mm.variance"] = flex.vec3_double(var_x, var_y, var_z)

    # use a ReflectionManager to exclude reflections too close to the spindle

    refman = LaueReflectionManager(obs_refs, experiments, outlier_detector=None)
    refman.finalise()

    # Redefine the reflection predictor to use the type expected by the Target class
    ref_predictor = LaueExperimentsPredictor(experiments)

    # keep only those reflections that pass inclusion criteria and have predictions
    reflections = refman.get_matches()

    # get analytical gradients
    an_grads = pred_param.get_gradients(reflections)

    # get finite difference gradients
    p_vals = pred_param.get_param_vals()
    deltas = [1.0e-7] * len(p_vals)

    for i, delta in enumerate(deltas):
        val = p_vals[i]

        p_vals[i] -= delta / 2.0
        pred_param.set_param_vals(p_vals)

        ref_predictor(reflections)

        rev_state = reflections["xyzcal.mm"].deep_copy()
        rev_wavelengths = reflections["wavelength_cal"].deep_copy()

        p_vals[i] += delta
        pred_param.set_param_vals(p_vals)

        ref_predictor(reflections)

        fwd_state = reflections["xyzcal.mm"].deep_copy()
        fwd_wavelengths = reflections["wavelength_cal"].deep_copy()
        p_vals[i] = val

        fd = fwd_state - rev_state
        wavelength_grads = fwd_wavelengths - rev_wavelengths
        x_grads, y_grads, _ = fd.parts()
        x_grads /= delta
        y_grads /= delta
        wavelength_grads /= delta

        # compare with analytical calculation
        assert x_grads == pytest.approx(an_grads[i]["dX_dp"], abs=5.0e-5)
        assert y_grads == pytest.approx(an_grads[i]["dY_dp"], abs=5.0e-5)
        assert wavelength_grads == pytest.approx(
            an_grads[i]["dwavelength_dp"], abs=5.0e-7
        )

    # return to the initial state
    pred_param.set_param_vals(p_vals)
