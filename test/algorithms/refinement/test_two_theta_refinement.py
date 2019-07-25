#
#  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""
Test refinement of a crystal unit cell using a two theta target.

"""

from __future__ import absolute_import, division, print_function
import os
from libtbx.test_utils import approx_equal
from math import pi
from copy import deepcopy

from dxtbx.model.experiment_list import ExperimentList, Experiment
from dials.algorithms.refinement.two_theta_refiner import (
    TwoThetaTarget,
    TwoThetaReflectionManager,
    TwoThetaExperimentsPredictor,
    TwoThetaPredictionParameterisation,
)


def generate_reflections(experiments):

    from dials.algorithms.spot_prediction import IndexGenerator
    from dials.algorithms.refinement.prediction.managed_predictors import (
        ScansRayPredictor,
        ScansExperimentsPredictor,
    )
    from dials.algorithms.spot_prediction import ray_intersection
    from cctbx.sgtbx import space_group, space_group_symbols
    from scitbx.array_family import flex

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

    # Invent some variances for the centroid positions of the simulated data
    im_width = 0.1 * pi / 180.0
    px_size = detector[0].get_pixel_size()
    var_x = flex.double(len(obs_refs), (px_size[0] / 2.0) ** 2)
    var_y = flex.double(len(obs_refs), (px_size[1] / 2.0) ** 2)
    var_phi = flex.double(len(obs_refs), (im_width / 2.0) ** 2)
    obs_refs["xyzobs.mm.variance"] = flex.vec3_double(var_x, var_y, var_phi)

    return obs_refs, ref_predictor


def test_fd_derivatives():
    """Test derivatives of the prediction equation"""

    from libtbx.phil import parse

    # Import model builder
    from dials.test.algorithms.refinement.setup_geometry import Extract

    # Create models
    overrides = """geometry.parameters.crystal.a.length.range = 10 50
  geometry.parameters.crystal.b.length.range = 10 50
  geometry.parameters.crystal.c.length.range = 10 50"""
    master_phil = parse(
        """
      include scope dials.test.algorithms.refinement.geometry_phil
      """,
        process_includes=True,
    )
    models = Extract(master_phil, overrides)

    mydetector = models.detector
    mygonio = models.goniometer
    mycrystal = models.crystal
    mybeam = models.beam

    # Build a mock scan for a 72 degree sweep
    from dxtbx.model import ScanFactory

    sf = ScanFactory()
    myscan = sf.make_scan(
        image_range=(1, 720),
        exposure_times=0.1,
        oscillation=(0, 0.1),
        epochs=list(range(720)),
        deg=True,
    )

    # Create a parameterisation of the crystal unit cell
    from dials.algorithms.refinement.parameterisation.crystal_parameters import (
        CrystalUnitCellParameterisation,
    )

    xluc_param = CrystalUnitCellParameterisation(mycrystal)

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

    # Build a prediction parameterisation for two theta prediction
    pred_param = TwoThetaPredictionParameterisation(
        experiments,
        detector_parameterisations=None,
        beam_parameterisations=None,
        xl_orientation_parameterisations=None,
        xl_unit_cell_parameterisations=[xluc_param],
    )

    # Generate some reflections
    obs_refs, ref_predictor = generate_reflections(experiments)

    # Build a ReflectionManager with overloads for handling 2theta residuals
    refman = TwoThetaReflectionManager(obs_refs, experiments, outlier_detector=None)

    # Build a TwoThetaExperimentsPredictor
    ref_predictor = TwoThetaExperimentsPredictor(experiments)

    # Make a target for the least squares 2theta residual
    target = TwoThetaTarget(experiments, ref_predictor, refman, pred_param)

    # Keep only reflections that pass inclusion criteria and have predictions
    reflections = refman.get_matches()

    # Get analytical gradients
    an_grads = pred_param.get_gradients(reflections)

    # Get finite difference gradients
    p_vals = pred_param.get_param_vals()
    deltas = [1.0e-7] * len(p_vals)

    for i in range(len(deltas)):

        val = p_vals[i]

        p_vals[i] -= deltas[i] / 2.0
        pred_param.set_param_vals(p_vals)

        target.predict()
        reflections = refman.get_matches()

        rev_state = reflections["2theta_resid"].deep_copy()

        p_vals[i] += deltas[i]
        pred_param.set_param_vals(p_vals)

        target.predict()
        reflections = refman.get_matches()

        fwd_state = reflections["2theta_resid"].deep_copy()
        p_vals[i] = val

        fd = fwd_state - rev_state
        fd /= deltas[i]

        # compare with analytical calculation
        assert approx_equal(fd, an_grads[i]["d2theta_dp"], eps=1.0e-6)

    # return to the initial state
    pred_param.set_param_vals(p_vals)


def test_refinement(dials_regression):
    """Test a refinement run"""

    # Get a beam and detector from a experiments. This one has a CS-PAD, but that
    # is irrelevant
    data_dir = os.path.join(dials_regression, "refinement_test_data", "hierarchy_test")
    experiments_path = os.path.join(data_dir, "datablock.json")
    assert os.path.exists(experiments_path)

    # load models
    from dxtbx.model.experiment_list import ExperimentListFactory

    experiments = ExperimentListFactory.from_serialized_format(
        experiments_path, check_format=False
    )
    im_set = experiments.imagesets()[0]
    detector = deepcopy(im_set.get_detector())
    beam = im_set.get_beam()

    # Invent a crystal, goniometer and scan for this test
    from dxtbx.model import Crystal

    crystal = Crystal(
        (40.0, 0.0, 0.0), (0.0, 40.0, 0.0), (0.0, 0.0, 40.0), space_group_symbol="P1"
    )
    orig_xl = deepcopy(crystal)

    from dxtbx.model import GoniometerFactory

    goniometer = GoniometerFactory.known_axis((1.0, 0.0, 0.0))

    # Build a mock scan for a 180 degree sweep
    from dxtbx.model import ScanFactory

    sf = ScanFactory()
    scan = sf.make_scan(
        image_range=(1, 1800),
        exposure_times=0.1,
        oscillation=(0, 0.1),
        epochs=list(range(1800)),
        deg=True,
    )
    sweep_range = scan.get_oscillation_range(deg=False)
    im_width = scan.get_oscillation(deg=False)[1]
    assert sweep_range == (0.0, pi)
    assert approx_equal(im_width, 0.1 * pi / 180.0)

    # Build an experiment list
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

    # simulate some reflections
    refs, _ = generate_reflections(experiments)

    # change unit cell a bit (=0.1 Angstrom length upsets, 0.1 degree of
    # alpha and beta angles)
    from dials.algorithms.refinement.parameterisation.crystal_parameters import (
        CrystalUnitCellParameterisation,
    )

    xluc_param = CrystalUnitCellParameterisation(crystal)
    cell_params = crystal.get_unit_cell().parameters()
    cell_params = [a + b for a, b in zip(cell_params, [0.1, -0.1, 0.1, 0.1, -0.1, 0.0])]
    from cctbx.uctbx import unit_cell
    from rstbx.symmetry.constraints.parameter_reduction import symmetrize_reduce_enlarge
    from scitbx import matrix

    new_uc = unit_cell(cell_params)
    newB = matrix.sqr(new_uc.fractionalization_matrix()).transpose()
    S = symmetrize_reduce_enlarge(crystal.get_space_group())
    S.set_orientation(orientation=newB)
    X = tuple([e * 1.0e5 for e in S.forward_independent_parameters()])
    xluc_param.set_param_vals(X)

    # reparameterise the crystal at the perturbed geometry
    xluc_param = CrystalUnitCellParameterisation(crystal)

    # Dummy parameterisations for other models
    beam_param = None
    xlo_param = None
    det_param = None

    # parameterisation of the prediction equation
    from dials.algorithms.refinement.parameterisation.parameter_report import (
        ParameterReporter,
    )

    pred_param = TwoThetaPredictionParameterisation(
        experiments, det_param, beam_param, xlo_param, [xluc_param]
    )
    param_reporter = ParameterReporter(det_param, beam_param, xlo_param, [xluc_param])

    # reflection manager
    refman = TwoThetaReflectionManager(refs, experiments, nref_per_degree=20)

    # reflection predictor
    ref_predictor = TwoThetaExperimentsPredictor(experiments)

    # target function
    target = TwoThetaTarget(experiments, ref_predictor, refman, pred_param)

    # minimisation engine
    from dials.algorithms.refinement.engine import (
        LevenbergMarquardtIterations as Refinery,
    )

    refinery = Refinery(
        target=target,
        prediction_parameterisation=pred_param,
        log=None,
        max_iterations=20,
    )

    # Refiner
    from dials.algorithms.refinement.refiner import Refiner

    refiner = Refiner(
        experiments=experiments,
        pred_param=pred_param,
        param_reporter=param_reporter,
        refman=refman,
        target=target,
        refinery=refinery,
    )
    refiner.run()

    # compare crystal with original crystal
    refined_xl = refiner.get_experiments()[0].crystal

    # print refined_xl
    assert refined_xl.is_similar_to(
        orig_xl, uc_rel_length_tolerance=0.001, uc_abs_angle_tolerance=0.01
    )

    # print "Unit cell esds:"
    # print refined_xl.get_cell_parameter_sd()
