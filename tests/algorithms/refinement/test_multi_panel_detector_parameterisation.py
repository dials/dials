"""
Test refinement of beam, detector and crystal orientation parameters using
generated reflection positions from ideal geometry, repeating tests with both a
single panel detector, and a geometrically identical 3x3 panel detector,
ensuring the results are the same.
"""


from __future__ import annotations

from collections import namedtuple
from math import pi

import pytest

from cctbx.sgtbx import space_group, space_group_symbols
from cctbx.uctbx import unit_cell
from dxtbx.model import Detector, Panel, ScanFactory
from dxtbx.model.experiment_list import Experiment, ExperimentList
from libtbx.phil import parse
from libtbx.test_utils import approx_equal
from rstbx.symmetry.constraints.parameter_reduction import symmetrize_reduce_enlarge
from scitbx import matrix

from dials.algorithms.refinement.parameterisation.beam_parameters import (
    BeamParameterisation,
)
from dials.algorithms.refinement.parameterisation.crystal_parameters import (
    CrystalOrientationParameterisation,
    CrystalUnitCellParameterisation,
)
from dials.algorithms.refinement.parameterisation.detector_parameters import (
    DetectorParameterisationMultiPanel,
    DetectorParameterisationSinglePanel,
    PyDetectorParameterisationMultiPanel,
)
from dials.algorithms.refinement.parameterisation.prediction_parameters import (
    XYPhiPredictionParameterisation,
)
from dials.algorithms.refinement.prediction.managed_predictors import (
    ScansExperimentsPredictor,
    ScansRayPredictor,
)
from dials.algorithms.refinement.reflection_manager import ReflectionManager
from dials.algorithms.refinement.target import (
    LeastSquaresPositionalResidualWithRmsdCutoff,
)
from dials.algorithms.spot_prediction import IndexGenerator, ray_intersection
from dials.array_family import flex

from . import geometry_phil, minimiser_phil, setup_geometry, setup_minimiser


def make_panel_in_array(array_elt, reference_panel):
    """Helper function to make a panel in a coplanar array with each panel size
    1/3 that of a reference panel"""

    px_size = tuple((e / 3.0) for e in reference_panel.get_pixel_size())
    ref_panel_size = reference_panel.get_image_size_mm()
    x_shift = array_elt[0] * ref_panel_size[0] / 3.0
    y_shift = array_elt[1] * ref_panel_size[1] / 3.0
    origin = (
        matrix.col(reference_panel.get_origin())
        + x_shift * matrix.col(reference_panel.get_fast_axis())
        + y_shift * matrix.col(reference_panel.get_slow_axis())
    )
    return Panel(
        type="PAD",
        name="Panel",
        fast_axis=reference_panel.get_fast_axis(),
        slow_axis=reference_panel.get_slow_axis(),
        origin=origin,
        pixel_size=px_size,
        image_size=reference_panel.get_image_size(),
        trusted_range=(0, 1.0e6),
        thickness=0.0,
        material="",
    )


# Setup experimental models
master_phil = parse(
    f"""
    {geometry_phil}
    {minimiser_phil}
    """
)


@pytest.fixture(scope="session")
def init_test():

    models = setup_geometry.Extract(master_phil)

    single_panel_detector = models.detector
    gonio = models.goniometer
    crystal = models.crystal
    beam = models.beam

    # Make a 3x3 multi panel detector filling the same space as the existing
    # single panel detector. Each panel of the multi-panel detector has pixels
    # with 1/3 the length dimensions of the single panel.
    multi_panel_detector = Detector()
    for x in range(3):
        for y in range(3):
            new_panel = make_panel_in_array((x, y), single_panel_detector[0])
            multi_panel_detector.add_panel(new_panel)

    # Build a mock scan for a 180 degree sequence
    sf = ScanFactory()
    scan = sf.make_scan(
        image_range=(1, 1800),
        exposure_times=0.1,
        oscillation=(0, 0.1),
        epochs=list(range(1800)),
        deg=True,
    )
    sequence_range = scan.get_oscillation_range(deg=False)
    im_width = scan.get_oscillation(deg=False)[1]
    assert sequence_range == (0.0, pi)
    assert approx_equal(im_width, 0.1 * pi / 180.0)

    # Build ExperimentLists
    experiments_single_panel = ExperimentList()
    experiments_multi_panel = ExperimentList()
    experiments_single_panel.append(
        Experiment(
            beam=beam,
            detector=single_panel_detector,
            goniometer=gonio,
            scan=scan,
            crystal=crystal,
            imageset=None,
        )
    )
    experiments_multi_panel.append(
        Experiment(
            beam=beam,
            detector=multi_panel_detector,
            goniometer=gonio,
            scan=scan,
            crystal=crystal,
            imageset=None,
        )
    )

    # Generate some reflections

    # All indices in a 2.0 Angstrom sphere
    resolution = 2.0
    index_generator = IndexGenerator(
        crystal.get_unit_cell(),
        space_group(space_group_symbols(1).hall()).type(),
        resolution,
    )
    indices = index_generator.to_array()

    # for the reflection predictor, it doesn't matter which experiment list is
    # passed, as the detector is not used
    ref_predictor = ScansRayPredictor(
        experiments_single_panel, scan.get_oscillation_range(deg=False)
    )

    # get two sets of identical reflections
    obs_refs_single = ref_predictor(indices)
    obs_refs_multi = ref_predictor(indices)
    for r1, r2 in zip(obs_refs_single.rows(), obs_refs_multi.rows()):
        assert r1["s1"] == r2["s1"]

    # get the panel intersections
    sel = ray_intersection(single_panel_detector, obs_refs_single)
    obs_refs_single = obs_refs_single.select(sel)
    sel = ray_intersection(multi_panel_detector, obs_refs_multi)
    obs_refs_multi = obs_refs_multi.select(sel)
    assert len(obs_refs_single) == len(obs_refs_multi)

    # Set 'observed' centroids from the predicted ones
    obs_refs_single["xyzobs.mm.value"] = obs_refs_single["xyzcal.mm"]
    obs_refs_multi["xyzobs.mm.value"] = obs_refs_multi["xyzcal.mm"]

    # Invent some variances for the centroid positions of the simulated data
    im_width = 0.1 * pi / 180.0
    px_size = single_panel_detector[0].get_pixel_size()
    var_x = flex.double(len(obs_refs_single), (px_size[0] / 2.0) ** 2)
    var_y = flex.double(len(obs_refs_single), (px_size[1] / 2.0) ** 2)
    var_phi = flex.double(len(obs_refs_single), (im_width / 2.0) ** 2)

    # set the variances and frame numbers
    obs_refs_single["xyzobs.mm.variance"] = flex.vec3_double(var_x, var_y, var_phi)
    obs_refs_multi["xyzobs.mm.variance"] = flex.vec3_double(var_x, var_y, var_phi)

    # Add in flags and ID columns by copying into standard reflection tables
    tmp = flex.reflection_table.empty_standard(len(obs_refs_single))
    tmp.update(obs_refs_single)
    obs_refs_single = tmp
    tmp = flex.reflection_table.empty_standard(len(obs_refs_multi))
    tmp.update(obs_refs_multi)
    obs_refs_multi = tmp

    test_data = namedtuple(
        "test_data",
        [
            "experiments_single_panel",
            "experiments_multi_panel",
            "observations_single_panel",
            "observations_multi_panel",
        ],
    )

    return test_data(
        experiments_single_panel,
        experiments_multi_panel,
        obs_refs_single,
        obs_refs_multi,
    )


def test(init_test):

    single_panel_detector = init_test.experiments_single_panel.detectors()[0]
    multi_panel_detector = init_test.experiments_multi_panel.detectors()[0]
    beam = init_test.experiments_single_panel.beams()[0]
    gonio = init_test.experiments_single_panel.goniometers()[0]
    crystal = init_test.experiments_single_panel.crystals()[0]

    # Parameterise the models
    det_param = DetectorParameterisationSinglePanel(single_panel_detector)
    s0_param = BeamParameterisation(beam, gonio)
    xlo_param = CrystalOrientationParameterisation(crystal)
    xluc_param = CrystalUnitCellParameterisation(crystal)

    multi_det_param = DetectorParameterisationMultiPanel(multi_panel_detector, beam)

    # Fix beam to the X-Z plane (imgCIF geometry), fix wavelength
    s0_param.set_fixed([True, False, True])

    # Link model parameterisations together into a parameterisation of the
    # prediction equation, first for the single panel detector
    pred_param = XYPhiPredictionParameterisation(
        init_test.experiments_single_panel,
        [det_param],
        [s0_param],
        [xlo_param],
        [xluc_param],
    )

    # ... and now for the multi-panel detector
    pred_param2 = XYPhiPredictionParameterisation(
        init_test.experiments_multi_panel,
        [multi_det_param],
        [s0_param],
        [xlo_param],
        [xluc_param],
    )

    ################################
    # Apply known parameter shifts #
    ################################

    # shift detectors by 1.0 mm each translation and 2 mrad each rotation
    det_p_vals = det_param.get_param_vals()
    p_vals = [a + b for a, b in zip(det_p_vals, [1.0, 1.0, 1.0, 2.0, 2.0, 2.0])]
    det_param.set_param_vals(p_vals)

    multi_det_p_vals = multi_det_param.get_param_vals()
    p_vals = [a + b for a, b in zip(multi_det_p_vals, [1.0, 1.0, 1.0, 2.0, 2.0, 2.0])]
    multi_det_param.set_param_vals(p_vals)

    # shift beam by 2 mrad in free axis
    s0_p_vals = s0_param.get_param_vals()
    p_vals = list(s0_p_vals)

    p_vals[0] += 2.0
    s0_param.set_param_vals(p_vals)

    # rotate crystal a bit (=2 mrad each rotation)
    xlo_p_vals = xlo_param.get_param_vals()
    p_vals = [a + b for a, b in zip(xlo_p_vals, [2.0, 2.0, 2.0])]
    xlo_param.set_param_vals(p_vals)

    # change unit cell a bit (=0.1 Angstrom length upsets, 0.1 degree of
    # gamma angle)
    xluc_p_vals = xluc_param.get_param_vals()
    cell_params = crystal.get_unit_cell().parameters()
    cell_params = [a + b for a, b in zip(cell_params, [0.1, 0.1, 0.1, 0.0, 0.0, 0.1])]
    new_uc = unit_cell(cell_params)
    newB = matrix.sqr(new_uc.fractionalization_matrix()).transpose()
    S = symmetrize_reduce_enlarge(crystal.get_space_group())
    S.set_orientation(orientation=newB)
    X = tuple([e * 1.0e5 for e in S.forward_independent_parameters()])
    xluc_param.set_param_vals(X)

    ###############################
    # Undo known parameter shifts #
    ###############################

    s0_param.set_param_vals(s0_p_vals)
    det_param.set_param_vals(det_p_vals)
    multi_det_param.set_param_vals(det_p_vals)
    xlo_param.set_param_vals(xlo_p_vals)
    xluc_param.set_param_vals(xluc_p_vals)

    #####################################
    # Select reflections for refinement #
    #####################################

    refman = ReflectionManager(
        init_test.observations_single_panel, init_test.experiments_single_panel
    )
    refman2 = ReflectionManager(
        init_test.observations_multi_panel, init_test.experiments_multi_panel
    )

    ###############################
    # Set up the target functions #
    ###############################

    target = LeastSquaresPositionalResidualWithRmsdCutoff(
        init_test.experiments_single_panel,
        ScansExperimentsPredictor(init_test.experiments_single_panel),
        refman,
        pred_param,
        restraints_parameterisation=None,
    )
    target2 = LeastSquaresPositionalResidualWithRmsdCutoff(
        init_test.experiments_multi_panel,
        ScansExperimentsPredictor(init_test.experiments_multi_panel),
        refman2,
        pred_param2,
        restraints_parameterisation=None,
    )

    #################################
    # Set up the refinement engines #
    #################################

    refiner = setup_minimiser.Extract(master_phil, target, pred_param).refiner
    refiner2 = setup_minimiser.Extract(master_phil, target2, pred_param2).refiner

    refiner.run()

    # reset parameters and run refinement with the multi panel detector
    s0_param.set_param_vals(s0_p_vals)
    multi_det_param.set_param_vals(det_p_vals)
    xlo_param.set_param_vals(xlo_p_vals)
    xluc_param.set_param_vals(xluc_p_vals)

    refiner2.run()

    # same number of steps
    assert refiner.get_num_steps() == refiner2.get_num_steps()

    # same rmsds
    for rmsd, rmsd2 in zip(refiner.history["rmsd"], refiner2.history["rmsd"]):
        assert approx_equal(rmsd, rmsd2)

    # same parameter values each step
    for params, params2 in zip(
        refiner.history["parameter_vector"], refiner.history["parameter_vector"]
    ):
        assert approx_equal(params, params2)


def test_equivalence_of_python_and_cpp_multipanel_algorithms(init_test):

    multi_panel_detector = init_test.experiments_multi_panel.detectors()[0]
    beam = init_test.experiments_single_panel.beams()[0]

    # Parameterise the models
    det_param1 = DetectorParameterisationMultiPanel(multi_panel_detector, beam)
    det_param2 = PyDetectorParameterisationMultiPanel(multi_panel_detector, beam)

    # shift detectors by 1.0 mm each translation and 2 mrad each rotation
    for dp in [det_param1, det_param2]:
        p_vals = dp.get_param_vals()
        p_vals = [a + b for a, b in zip(p_vals, [1.0, 1.0, 1.0, 2.0, 2.0, 2.0])]
        dp.set_param_vals(p_vals)
        dp.compose()

    for pnl in range(3):
        derivatives1 = det_param1.get_ds_dp(multi_state_elt=pnl)
        derivatives2 = det_param2.get_ds_dp(multi_state_elt=pnl)

        for a, b in zip(derivatives1, derivatives2):
            for i, j in zip(a, b):
                assert i == pytest.approx(j)
