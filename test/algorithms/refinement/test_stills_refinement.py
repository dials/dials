"""
A simple test of stills refinement using fake data.
Only the crystal is perturbed while the beam and detector are known.
"""

from __future__ import absolute_import, division, print_function


def test(args=[]):
    # Python and cctbx imports
    from math import pi
    from scitbx import matrix
    from libtbx.phil import parse
    from libtbx.test_utils import approx_equal

    # Import for surgery on reflection_tables
    from dials.array_family import flex

    # Get module to build models using PHIL
    import dials.test.algorithms.refinement.setup_geometry as setup_geometry

    # We will set up a mock scan and a mock experiment list
    from dxtbx.model import ScanFactory
    from dxtbx.model.experiment_list import ExperimentList, Experiment

    # Crystal parameterisations
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
        StillsExperimentsPredictor,
    )
    from dials.algorithms.spot_prediction import ray_intersection
    from cctbx.sgtbx import space_group, space_group_symbols

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

    # build models, with a larger crystal than default in order to get enough
    # reflections on the 'still' image
    param = """
  geometry.parameters.crystal.a.length.range=40 50;
  geometry.parameters.crystal.b.length.range=40 50;
  geometry.parameters.crystal.c.length.range=40 50;
  geometry.parameters.random_seed = 42"""
    models = setup_geometry.Extract(
        master_phil, cmdline_args=args, local_overrides=param
    )

    crystal = models.crystal
    mydetector = models.detector
    mygonio = models.goniometer
    mybeam = models.beam

    # Build a mock scan for a 1.5 degree wedge. Only used for generating indices near
    # the Ewald sphere
    sf = ScanFactory()
    myscan = sf.make_scan(
        image_range=(1, 1),
        exposure_times=0.1,
        oscillation=(0, 1.5),
        epochs=list(range(1)),
        deg=True,
    )
    sweep_range = myscan.get_oscillation_range(deg=False)
    im_width = myscan.get_oscillation(deg=False)[1]
    assert approx_equal(im_width, 1.5 * pi / 180.0)

    # Build experiment lists
    stills_experiments = ExperimentList()
    stills_experiments.append(
        Experiment(beam=mybeam, detector=mydetector, crystal=crystal, imageset=None)
    )
    scans_experiments = ExperimentList()
    scans_experiments.append(
        Experiment(
            beam=mybeam,
            detector=mydetector,
            crystal=crystal,
            goniometer=mygonio,
            scan=myscan,
            imageset=None,
        )
    )

    ##########################################################
    # Parameterise the models (only for perturbing geometry) #
    ##########################################################

    xlo_param = CrystalOrientationParameterisation(crystal)
    xluc_param = CrystalUnitCellParameterisation(crystal)

    ################################
    # Apply known parameter shifts #
    ################################

    # rotate crystal (=5 mrad each rotation)
    xlo_p_vals = []
    p_vals = xlo_param.get_param_vals()
    xlo_p_vals.append(p_vals)
    new_p_vals = [a + b for a, b in zip(p_vals, [5.0, 5.0, 5.0])]
    xlo_param.set_param_vals(new_p_vals)

    # change unit cell (=1.0 Angstrom length upsets, 0.5 degree of
    # gamma angle)
    xluc_p_vals = []
    p_vals = xluc_param.get_param_vals()
    xluc_p_vals.append(p_vals)
    cell_params = crystal.get_unit_cell().parameters()
    cell_params = [a + b for a, b in zip(cell_params, [1.0, 1.0, -1.0, 0.0, 0.0, 0.5])]
    new_uc = unit_cell(cell_params)
    newB = matrix.sqr(new_uc.fractionalization_matrix()).transpose()
    S = symmetrize_reduce_enlarge(crystal.get_space_group())
    S.set_orientation(orientation=newB)
    X = tuple([e * 1.0e5 for e in S.forward_independent_parameters()])
    xluc_param.set_param_vals(X)

    # keep track of the target crystal model to compare with refined
    from copy import deepcopy

    target_crystal = deepcopy(crystal)

    #############################
    # Generate some reflections #
    #############################

    # All indices in a 2.0 Angstrom sphere for crystal
    resolution = 2.0
    index_generator = IndexGenerator(
        crystal.get_unit_cell(),
        space_group(space_group_symbols(1).hall()).type(),
        resolution,
    )
    indices = index_generator.to_array()

    # Build a ray predictor and predict rays close to the Ewald sphere by using
    # the narrow rotation scan
    ref_predictor = ScansRayPredictor(scans_experiments, sweep_range)
    obs_refs = ref_predictor(indices, experiment_id=0)

    # Take only those rays that intersect the detector
    intersects = ray_intersection(mydetector, obs_refs)
    obs_refs = obs_refs.select(intersects)

    # Add in flags and ID columns by copying into standard reflection table
    tmp = flex.reflection_table.empty_standard(len(obs_refs))
    tmp.update(obs_refs)
    obs_refs = tmp

    # Invent some variances for the centroid positions of the simulated data
    im_width = 0.1 * pi / 180.0
    px_size = mydetector[0].get_pixel_size()
    var_x = flex.double(len(obs_refs), (px_size[0] / 2.0) ** 2)
    var_y = flex.double(len(obs_refs), (px_size[1] / 2.0) ** 2)
    var_phi = flex.double(len(obs_refs), (im_width / 2.0) ** 2)
    obs_refs["xyzobs.mm.variance"] = flex.vec3_double(var_x, var_y, var_phi)

    # Re-predict using the stills reflection predictor
    stills_ref_predictor = StillsExperimentsPredictor(stills_experiments)
    obs_refs_stills = stills_ref_predictor(obs_refs)

    # Set 'observed' centroids from the predicted ones
    obs_refs_stills["xyzobs.mm.value"] = obs_refs_stills["xyzcal.mm"]

    ###############################
    # Undo known parameter shifts #
    ###############################

    xlo_param.set_param_vals(xlo_p_vals[0])
    xluc_param.set_param_vals(xluc_p_vals[0])

    # make a refiner
    from dials.algorithms.refinement.refiner import phil_scope

    params = phil_scope.fetch(source=parse("")).extract()

    # Change this to get a plot
    do_plot = False
    if do_plot:
        params.refinement.refinery.journal.track_parameter_correlation = True

    from dials.algorithms.refinement.refiner import RefinerFactory

    # decrease bin_size_fraction to terminate on RMSD convergence
    params.refinement.target.bin_size_fraction = 0.01
    params.refinement.parameterisation.beam.fix = "all"
    params.refinement.parameterisation.detector.fix = "all"
    refiner = RefinerFactory.from_parameters_data_experiments(
        params, obs_refs_stills, stills_experiments
    )

    # run refinement
    history = refiner.run()

    # regression tests
    assert len(history["rmsd"]) == 9

    refined_crystal = refiner.get_experiments()[0].crystal
    uc1 = refined_crystal.get_unit_cell()
    uc2 = target_crystal.get_unit_cell()
    assert uc1.is_similar_to(uc2)

    if do_plot:
        plt = refiner.parameter_correlation_plot(
            len(history["parameter_correlation"]) - 1
        )
        plt.show()
