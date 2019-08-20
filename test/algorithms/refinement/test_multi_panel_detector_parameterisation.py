"""
Using code copied from tst_orientation_refinement.py, test refinement of beam,
detector and crystal orientation parameters using generated reflection positions
from ideal geometry, repeating tests with both a single panel detector, and a
geometrically identical 3x3 panel detector, ensuring the results are the same.

Control of the experimental model and choice of minimiser is done via
PHIL, which means we can do, for example:

cctbx.python tst_multi_panel_detector_parameterisation.py "random_seed=3"
"""

# Python and cctbx imports
from __future__ import absolute_import, division, print_function
from math import pi
from scitbx import matrix
from dials.array_family import flex
from libtbx.phil import parse
from libtbx.test_utils import approx_equal

# Get modules to build models and minimiser using PHIL
import dials.test.algorithms.refinement.setup_geometry as setup_geometry
import dials.test.algorithms.refinement.setup_minimiser as setup_minimiser

# Get the models to build the multi panel detector
from dxtbx.model import Panel, Detector

# We will set up a mock scan and a mock experiment list
from dxtbx.model import ScanFactory
from dxtbx.model.experiment_list import ExperimentList, Experiment

# Model parameterisations
from dials.algorithms.refinement.parameterisation.detector_parameters import (
    DetectorParameterisationSinglePanel,
    DetectorParameterisationMultiPanel,
)
from dials.algorithms.refinement.parameterisation.beam_parameters import (
    BeamParameterisation,
)
from dials.algorithms.refinement.parameterisation.crystal_parameters import (
    CrystalOrientationParameterisation,
    CrystalUnitCellParameterisation,
)

# Symmetry constrained parameterisation for the unit cell
from cctbx.uctbx import unit_cell
from rstbx.symmetry.constraints.parameter_reduction import symmetrize_reduce_enlarge

# Reflection prediction
from dials.algorithms.spot_prediction import IndexGenerator
from dials.algorithms.refinement.prediction.managed_predictors import (
    ScansRayPredictor,
    ScansExperimentsPredictor,
)
from dials.algorithms.spot_prediction import ray_intersection
from cctbx.sgtbx import space_group, space_group_symbols

# Parameterisation of the prediction equation
from dials.algorithms.refinement.parameterisation.prediction_parameters import (
    XYPhiPredictionParameterisation,
)

# Imports for the target function
from dials.algorithms.refinement.target import (
    LeastSquaresPositionalResidualWithRmsdCutoff,
)
from dials.algorithms.refinement.reflection_manager import ReflectionManager

###################
# Local functions #
###################


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


def test(args=[]):
    #############################
    # Setup experimental models #
    #############################
    master_phil = parse(
        """
      include scope dials.test.algorithms.refinement.geometry_phil
      include scope dials.test.algorithms.refinement.minimiser_phil
      """,
        process_includes=True,
    )

    models = setup_geometry.Extract(master_phil, cmdline_args=args)

    single_panel_detector = models.detector
    mygonio = models.goniometer
    mycrystal = models.crystal
    mybeam = models.beam

    # Make a 3x3 multi panel detector filling the same space as the existing
    # single panel detector. Each panel of the multi-panel detector has pixels with
    # 1/3 the length dimensions of the single panel.

    multi_panel_detector = Detector()
    for x in range(3):
        for y in range(3):
            new_panel = make_panel_in_array((x, y), single_panel_detector[0])
            multi_panel_detector.add_panel(new_panel)

    # Build a mock scan for a 180 degree sweep
    sf = ScanFactory()
    myscan = sf.make_scan(
        image_range=(1, 1800),
        exposure_times=0.1,
        oscillation=(0, 0.1),
        epochs=list(range(1800)),
        deg=True,
    )
    sweep_range = myscan.get_oscillation_range(deg=False)
    im_width = myscan.get_oscillation(deg=False)[1]
    assert sweep_range == (0.0, pi)
    assert approx_equal(im_width, 0.1 * pi / 180.0)

    # Build ExperimentLists
    experiments_single_panel = ExperimentList()
    experiments_multi_panel = ExperimentList()
    experiments_single_panel.append(
        Experiment(
            beam=mybeam,
            detector=single_panel_detector,
            goniometer=mygonio,
            scan=myscan,
            crystal=mycrystal,
            imageset=None,
        )
    )
    experiments_multi_panel.append(
        Experiment(
            beam=mybeam,
            detector=multi_panel_detector,
            goniometer=mygonio,
            scan=myscan,
            crystal=mycrystal,
            imageset=None,
        )
    )

    ###########################
    # Parameterise the models #
    ###########################

    det_param = DetectorParameterisationSinglePanel(single_panel_detector)
    s0_param = BeamParameterisation(mybeam, mygonio)
    xlo_param = CrystalOrientationParameterisation(mycrystal)
    xluc_param = CrystalUnitCellParameterisation(mycrystal)

    multi_det_param = DetectorParameterisationMultiPanel(multi_panel_detector, mybeam)

    # Fix beam to the X-Z plane (imgCIF geometry), fix wavelength
    s0_param.set_fixed([True, False, True])

    # Fix crystal parameters
    # xluc_param.set_fixed([True, True, True, True, True, True])

    ########################################################################
    # Link model parameterisations together into a parameterisation of the #
    # prediction equation                                                  #
    ########################################################################

    pred_param = XYPhiPredictionParameterisation(
        experiments_single_panel, [det_param], [s0_param], [xlo_param], [xluc_param]
    )

    pred_param2 = XYPhiPredictionParameterisation(
        experiments_multi_panel,
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
    cell_params = mycrystal.get_unit_cell().parameters()
    cell_params = [a + b for a, b in zip(cell_params, [0.1, 0.1, 0.1, 0.0, 0.0, 0.1])]
    new_uc = unit_cell(cell_params)
    newB = matrix.sqr(new_uc.fractionalization_matrix()).transpose()
    S = symmetrize_reduce_enlarge(mycrystal.get_space_group())
    S.set_orientation(orientation=newB)
    X = tuple([e * 1.0e5 for e in S.forward_independent_parameters()])
    xluc_param.set_param_vals(X)

    #############################
    # Generate some reflections #
    #############################

    # All indices in a 2.0 Angstrom sphere
    resolution = 2.0
    index_generator = IndexGenerator(
        mycrystal.get_unit_cell(),
        space_group(space_group_symbols(1).hall()).type(),
        resolution,
    )
    indices = index_generator.to_array()

    # for the reflection predictor, it doesn't matter which experiment list is
    # passed, as the detector is not used
    ref_predictor = ScansRayPredictor(experiments_single_panel, sweep_range)

    # get two sets of identical reflections
    obs_refs = ref_predictor(indices)
    obs_refs2 = ref_predictor(indices)
    for r1, r2 in zip(obs_refs, obs_refs2):
        assert r1["s1"] == r2["s1"]

    # get the panel intersections
    sel = ray_intersection(single_panel_detector, obs_refs)
    obs_refs = obs_refs.select(sel)
    sel = ray_intersection(multi_panel_detector, obs_refs2)
    obs_refs2 = obs_refs2.select(sel)
    assert len(obs_refs) == len(obs_refs2)

    # Set 'observed' centroids from the predicted ones
    obs_refs["xyzobs.mm.value"] = obs_refs["xyzcal.mm"]
    obs_refs2["xyzobs.mm.value"] = obs_refs2["xyzcal.mm"]

    # Invent some variances for the centroid positions of the simulated data
    im_width = 0.1 * pi / 180.0
    px_size = single_panel_detector[0].get_pixel_size()
    var_x = flex.double(len(obs_refs), (px_size[0] / 2.0) ** 2)
    var_y = flex.double(len(obs_refs), (px_size[1] / 2.0) ** 2)
    var_phi = flex.double(len(obs_refs), (im_width / 2.0) ** 2)

    # set the variances and frame numbers
    obs_refs["xyzobs.mm.variance"] = flex.vec3_double(var_x, var_y, var_phi)
    obs_refs2["xyzobs.mm.variance"] = flex.vec3_double(var_x, var_y, var_phi)

    # Add in flags and ID columns by copying into standard reflection tables
    tmp = flex.reflection_table.empty_standard(len(obs_refs))
    tmp.update(obs_refs)
    obs_refs = tmp
    tmp = flex.reflection_table.empty_standard(len(obs_refs2))
    tmp.update(obs_refs2)
    obs_refs2 = tmp

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

    refman = ReflectionManager(obs_refs, experiments_single_panel)
    refman2 = ReflectionManager(obs_refs, experiments_multi_panel)

    ###############################
    # Set up the target functions #
    ###############################

    mytarget = LeastSquaresPositionalResidualWithRmsdCutoff(
        experiments_single_panel,
        ScansExperimentsPredictor(experiments_single_panel),
        refman,
        pred_param,
        restraints_parameterisation=None,
    )
    mytarget2 = LeastSquaresPositionalResidualWithRmsdCutoff(
        experiments_multi_panel,
        ScansExperimentsPredictor(experiments_multi_panel),
        refman2,
        pred_param2,
        restraints_parameterisation=None,
    )

    #################################
    # Set up the refinement engines #
    #################################

    refiner = setup_minimiser.Extract(
        master_phil, mytarget, pred_param, cmdline_args=args
    ).refiner
    refiner2 = setup_minimiser.Extract(
        master_phil, mytarget2, pred_param2, cmdline_args=args
    ).refiner

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
