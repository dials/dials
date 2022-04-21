"""Test analytical calculation of gradients of the target function versus finite
difference calculations"""


from __future__ import annotations


def test(args=[]):
    # Python and cctbx imports
    import random
    from math import pi

    from cctbx.sgtbx import space_group, space_group_symbols

    # We will set up a mock scan and a mock experiment list
    from dxtbx.model import ScanFactory
    from dxtbx.model.experiment_list import Experiment, ExperimentList
    from libtbx.phil import parse
    from libtbx.test_utils import approx_equal
    from scitbx import matrix
    from scitbx.array_family import flex

    from dials.algorithms.refinement.parameterisation.beam_parameters import (
        BeamParameterisation,
    )
    from dials.algorithms.refinement.parameterisation.crystal_parameters import (
        CrystalOrientationParameterisation,
        CrystalUnitCellParameterisation,
    )

    # Model parameterisations
    from dials.algorithms.refinement.parameterisation.detector_parameters import (
        DetectorParameterisationSinglePanel,
    )

    # Parameterisation of the prediction equation
    from dials.algorithms.refinement.parameterisation.prediction_parameters import (
        XYPhiPredictionParameterisation,
    )
    from dials.algorithms.refinement.prediction.managed_predictors import (
        ScansExperimentsPredictor,
        ScansRayPredictor,
    )
    from dials.algorithms.refinement.reflection_manager import ReflectionManager

    # Imports for the target function
    from dials.algorithms.refinement.target import (
        LeastSquaresPositionalResidualWithRmsdCutoff,
    )

    # Reflection prediction
    from dials.algorithms.spot_prediction import IndexGenerator, ray_intersection

    from . import geometry_phil

    # Experimental model builder
    from .setup_geometry import Extract

    # Local functions
    def random_direction_close_to(vector, sd=0.5):
        return vector.rotate_around_origin(
            matrix.col((random.random(), random.random(), random.random())).normalize(),
            random.gauss(0, sd),
            deg=True,
        )

    #############################
    # Setup experimental models #
    #############################

    # make a small cell to speed up calculations
    overrides = """geometry.parameters.crystal.a.length.range = 10 15
  geometry.parameters.crystal.b.length.range = 10 15
  geometry.parameters.crystal.c.length.range = 10 15"""

    master_phil = parse(geometry_phil)

    models = Extract(master_phil, overrides, cmdline_args=args)

    mydetector = models.detector
    mygonio = models.goniometer
    mycrystal = models.crystal
    mybeam = models.beam

    # Build a mock scan for a 180 degree sequence of 0.1 degree images
    sf = ScanFactory()
    myscan = sf.make_scan(
        image_range=(1, 1800),
        exposure_times=0.1,
        oscillation=(0, 0.1),
        epochs=list(range(1800)),
        deg=True,
    )
    sequence_range = myscan.get_oscillation_range(deg=False)
    im_width = myscan.get_oscillation(deg=False)[1]
    assert sequence_range == (0.0, pi)
    assert approx_equal(im_width, 0.1 * pi / 180.0)

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

    ###########################
    # Parameterise the models #
    ###########################

    det_param = DetectorParameterisationSinglePanel(mydetector)
    s0_param = BeamParameterisation(mybeam, mygonio)
    xlo_param = CrystalOrientationParameterisation(mycrystal)
    xluc_param = CrystalUnitCellParameterisation(mycrystal)

    ########################################################################
    # Link model parameterisations together into a parameterisation of the #
    # prediction equation                                                  #
    ########################################################################

    pred_param = XYPhiPredictionParameterisation(
        experiments, [det_param], [s0_param], [xlo_param], [xluc_param]
    )

    ################################
    # Apply known parameter shifts #
    ################################

    # shift detector by 0.2 mm each translation and 2 mrad each rotation
    det_p_vals = det_param.get_param_vals()
    p_vals = [a + b for a, b in zip(det_p_vals, [2.0, 2.0, 2.0, 2.0, 2.0, 2.0])]
    det_param.set_param_vals(p_vals)

    # shift beam by 2 mrad in one axis
    s0_p_vals = s0_param.get_param_vals()
    p_vals = list(s0_p_vals)
    p_vals[1] += 2.0
    s0_param.set_param_vals(p_vals)

    # rotate crystal a bit (=2 mrad each rotation)
    xlo_p_vals = xlo_param.get_param_vals()
    p_vals = [a + b for a, b in zip(xlo_p_vals, [2.0, 2.0, 2.0])]
    xlo_param.set_param_vals(p_vals)

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
    im_width = 0.1 * pi / 180.0
    px_size = mydetector[0].get_pixel_size()
    var_x = flex.double(len(obs_refs), (px_size[0] / 2.0) ** 2)
    var_y = flex.double(len(obs_refs), (px_size[1] / 2.0) ** 2)
    var_phi = flex.double(len(obs_refs), (im_width / 2.0) ** 2)
    obs_refs["xyzobs.mm.variance"] = flex.vec3_double(var_x, var_y, var_phi)

    ###############################
    # Undo known parameter shifts #
    ###############################

    s0_param.set_param_vals(s0_p_vals)
    det_param.set_param_vals(det_p_vals)
    xlo_param.set_param_vals(xlo_p_vals)

    #####################################
    # Select reflections for refinement #
    #####################################

    refman = ReflectionManager(obs_refs, experiments)

    ##############################
    # Set up the target function #
    ##############################

    # Redefine the reflection predictor to use the type expected by the Target class
    ref_predictor = ScansExperimentsPredictor(experiments)

    mytarget = LeastSquaresPositionalResidualWithRmsdCutoff(
        experiments, ref_predictor, refman, pred_param, restraints_parameterisation=None
    )

    # get the functional and gradients
    mytarget.predict()
    L, dL_dp, curvs = mytarget.compute_functional_gradients_and_curvatures()

    ####################################
    # Do FD calculation for comparison #
    ####################################

    # function for calculating finite difference gradients of the target function
    def get_fd_gradients(target, pred_param, deltas):
        """Calculate centered finite difference gradients for each of the
        parameters of the target function.

        "deltas" must be a sequence of the same length as the parameter list, and
        contains the step size for the difference calculations for each parameter.
        """

        p_vals = pred_param.get_param_vals()
        assert len(deltas) == len(p_vals)
        fd_grad = []
        fd_curvs = []

        for i in range(len(deltas)):
            val = p_vals[i]

            p_vals[i] -= deltas[i] / 2.0
            pred_param.set_param_vals(p_vals)
            target.predict()

            rev_state = target.compute_functional_gradients_and_curvatures()

            p_vals[i] += deltas[i]
            pred_param.set_param_vals(p_vals)

            target.predict()

            fwd_state = target.compute_functional_gradients_and_curvatures()

            # finite difference estimation of first derivatives
            fd_grad.append((fwd_state[0] - rev_state[0]) / deltas[i])

            # finite difference estimation of curvatures, using the analytical
            # first derivatives
            fd_curvs.append((fwd_state[1][i] - rev_state[1][i]) / deltas[i])

            # set parameter back to centred value
            p_vals[i] = val

        # return to the initial state
        pred_param.set_param_vals(p_vals)

        return fd_grad, fd_curvs

    # test normalised differences between FD and analytical calculations
    fdgrads = get_fd_gradients(mytarget, pred_param, [1.0e-7] * len(pred_param))
    diffs = [a - b for a, b in zip(dL_dp, fdgrads[0])]
    norm_diffs = tuple([a / b for a, b in zip(diffs, fdgrads[0])])
    for e in norm_diffs:
        assert abs(e) < 0.001  # check differences less than 0.1%

    # test normalised differences between FD curvatures and analytical least
    # squares approximation. We don't expect this to be especially close
    if curvs:
        diffs = [a - b for a, b in zip(curvs, fdgrads[1])]
        norm_diffs = tuple([a / b for a, b in zip(diffs, fdgrads[1])])
        for e in norm_diffs:
            assert abs(e) < 0.1  # check differences less than 10%
