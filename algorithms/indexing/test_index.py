from __future__ import absolute_import, division, print_function

import collections
import glob
import os
import pytest

import procrunner
from cctbx import uctbx
from dxtbx.serialize import load
from dxtbx.model import ExperimentList
from dials.array_family import flex


def unit_cells_are_similar(
    uc1, uc2, relative_length_tolerance=0.01, absolute_angle_tolerance=1
):
    # see also uctbx.cpp unit_cell::is_similar_to()
    l1 = uc1.parameters()
    l2 = uc2.parameters()
    for i in range(3):
        if abs(min(l1[i], l2[i]) / max(l1[i], l2[i]) - 1) > relative_length_tolerance:
            return False
    for i in range(3, 6):
        if abs(l1[i] - l2[i]) > absolute_angle_tolerance:
            if abs(l1[i] - (180 - l2[i])) > absolute_angle_tolerance:
                return False
    return True


_indexing_result = collections.namedtuple(
    "indexing", ["indexed_reflections", "experiments", "rmsds"]
)


def run_indexing(
    reflections,
    experiment,
    working_directory,
    extra_args,
    expected_unit_cell,
    expected_rmsds,
    expected_hall_symbol,
    n_expected_lattices=1,
    relative_length_tolerance=0.005,
    absolute_angle_tolerance=0.5,
):
    commands = ["dials.index"]
    if isinstance(reflections, list):
        commands.extend(reflections)
    else:
        commands.append(reflections)
    if isinstance(experiment, list):
        commands.extend(experiment)
    else:
        commands.append(experiment)
    commands.extend(extra_args)

    result = procrunner.run(commands, working_directory=working_directory)
    assert not result.returncode and not result.stderr

    out_expts = working_directory.join("indexed.expt")
    out_refls = working_directory.join("indexed.refl")
    assert out_expts.check()
    assert out_refls.check()

    experiments_list = load.experiment_list(out_expts.strpath, check_format=False)
    assert len(experiments_list.crystals()) == n_expected_lattices
    indexed_reflections = flex.reflection_table.from_file(out_refls.strpath)
    indexed_reflections.assert_experiment_identifiers_are_consistent(experiments_list)
    rmsds = None

    for i, experiment in enumerate(experiments_list):
        assert unit_cells_are_similar(
            experiment.crystal.get_unit_cell(),
            expected_unit_cell,
            relative_length_tolerance=relative_length_tolerance,
            absolute_angle_tolerance=absolute_angle_tolerance,
        ), (
            experiment.crystal.get_unit_cell().parameters(),
            expected_unit_cell.parameters(),
        )
        sg = experiment.crystal.get_space_group()
        assert sg.type().hall_symbol() == expected_hall_symbol, (
            sg.type().hall_symbol(),
            expected_hall_symbol,
        )
        reflections = indexed_reflections.select(indexed_reflections["id"] == i)
        mi = reflections["miller_index"]
        assert (mi != (0, 0, 0)).count(False) == 0
        reflections = reflections.select(mi != (0, 0, 0))
        reflections = reflections.select(
            reflections.get_flags(reflections.flags.used_in_refinement)
        )
        assert len(reflections) > 0
        obs_x, obs_y, obs_z = reflections["xyzobs.mm.value"].parts()
        calc_x, calc_y, calc_z = reflections["xyzcal.mm"].parts()
        rmsd_x = flex.mean(flex.pow2(obs_x - calc_x)) ** 0.5
        rmsd_y = flex.mean(flex.pow2(obs_y - calc_y)) ** 0.5
        rmsd_z = flex.mean(flex.pow2(obs_z - calc_z)) ** 0.5
        rmsds = (rmsd_x, rmsd_y, rmsd_z)
        for actual, expected in zip(rmsds, expected_rmsds):
            assert actual <= expected, "%s %s" % (rmsds, expected_rmsds)
        assert experiment.identifier != ""
        expt = ExperimentList()
        expt.append(experiment)
        reflections.assert_experiment_identifiers_are_consistent(expt)

    return _indexing_result(indexed_reflections, experiments_list, rmsds)


def test_index_i04_weak_data_fft3d(dials_regression, tmpdir):
    # thaumatin
    data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")
    pickle_path = os.path.join(data_dir, "full.pickle")
    sequence_path = os.path.join(data_dir, "experiments_import.json")
    extra_args = [
        "bin_size_fraction=0.25",
        "image_range=1,20",
        "image_range=250,270",
        "image_range=520,540",
    ]
    expected_unit_cell = uctbx.unit_cell((57.7, 57.7, 149.8, 90, 90, 90))
    expected_rmsds = (0.05, 0.04, 0.0005)
    expected_hall_symbol = " P 1"

    run_indexing(
        pickle_path,
        sequence_path,
        tmpdir,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )


def test_index_trypsin_four_lattice_P212121(dials_regression, tmpdir):
    # synthetic trypsin multi-lattice dataset (4 lattices)
    data_dir = os.path.join(dials_regression, "indexing_test_data", "trypsin")
    pickle_path = os.path.join(data_dir, "P1_X6_1_2_3_4.pickle")
    sequence_path = os.path.join(data_dir, "experiments_P1_X6_1_2_3_4.json")
    extra_args = [
        "indexing.method=real_space_grid_search",
        "reflections_per_degree=10",
        "n_macro_cycles=5",
        "known_symmetry.unit_cell=54.3,58.3,66.5,90,90,90",
        "known_symmetry.space_group=P212121",
        "image_range=0,10",
        "beam.fix=all",
        "detector.fix=all",
        "max_cell=70",
    ]
    expected_unit_cell = uctbx.unit_cell((54.3, 58.3, 66.5, 90, 90, 90))
    expected_rmsds = (0.28, 0.30, 0.006)
    expected_hall_symbol = " P 2ac 2ab"
    n_expected_lattices = 1

    run_indexing(
        pickle_path,
        sequence_path,
        tmpdir,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
        n_expected_lattices=n_expected_lattices,
        relative_length_tolerance=0.02,
        absolute_angle_tolerance=1,
    )


def test_index_i04_weak_data_fft1d(dials_regression, tmpdir):
    # thaumatin
    data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")
    pickle_path = os.path.join(data_dir, "full.pickle")
    sequence_path = os.path.join(data_dir, "experiments_import.json")
    extra_args = [
        "n_macro_cycles=2",
        "indexing.method=fft1d",
        "bin_size_fraction=0.25",
        "image_range=1,20",
        "image_range=250,270",
        "image_range=520,540",
    ]
    expected_unit_cell = uctbx.unit_cell((57.7, 57.7, 149.9, 90, 90, 90))
    expected_rmsds = (0.06, 0.05, 0.0005)
    expected_hall_symbol = " P 1"

    run_indexing(
        pickle_path,
        sequence_path,
        tmpdir,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )


def test_index_trypsin_index_assignment_local(dials_regression, tmpdir):
    # synthetic trypsin multi-lattice dataset (3 lattices)
    data_dir = os.path.join(dials_regression, "indexing_test_data", "trypsin")
    pickle_path = os.path.join(data_dir, "P1_X6_1_2_3.pickle")
    sequence_path = os.path.join(data_dir, "experiments_P1_X6_1_2_3.json")
    extra_args = [
        "indexing.method=real_space_grid_search",
        "d_min_start=3",
        "n_macro_cycles=3",
        "known_symmetry.unit_cell=54.3,58.3,66.5,90,90,90",
        "known_symmetry.space_group=P212121",
        "image_range=0,10",
        "beam.fix=all",
        "detector.fix=all",
        "max_lattices=3",
        "index_assignment.method=local",
        "nearest_neighbours=50",
    ]

    expected_unit_cell = uctbx.unit_cell((54.3, 58.3, 66.5, 90, 90, 90))
    expected_rmsds = (0.33, 0.40, 0.0024)
    expected_hall_symbol = " P 2ac 2ab"
    n_expected_lattices = 3

    run_indexing(
        pickle_path,
        sequence_path,
        tmpdir,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
        n_expected_lattices=n_expected_lattices,
        relative_length_tolerance=0.02,
        absolute_angle_tolerance=1,
    )


def test_index_peak_search_clean(dials_regression, tmpdir):
    # test indexing from single image of i04_weak_data
    data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")
    pickle_path = os.path.join(data_dir, "first_image.pickle")
    sequence_path = os.path.join(data_dir, "experiments_import.json")
    extra_args = [
        "indexing.method=fft3d",
        "known_symmetry.space_group=P4",
        "known_symmetry.unit_cell=57.8,57.8,150,90,90,90",
        "peak_search=clean",
        "min_samples=15",
        "n_macro_cycles=4",
        "reciprocal_space_grid.d_min=4",
    ]

    expected_unit_cell = uctbx.unit_cell((57.8, 57.8, 150, 90, 90, 90))
    expected_rmsds = (0.06, 0.07, 0.003)
    expected_hall_symbol = " P 4"

    run_indexing(
        pickle_path,
        sequence_path,
        tmpdir,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )


@pytest.mark.parametrize("specify_unit_cell", [False, True])
def test_index_imosflm_tutorial(dials_regression, tmpdir, specify_unit_cell):
    # test on spots derived from imosflm tutorial data:
    # http://www.ccp4.ac.uk/courses/BCA2005/tutorials/dataproc-tutorial.html
    data_dir = os.path.join(dials_regression, "indexing_test_data", "imosflm_hg_mar")
    pickle_path = os.path.join(data_dir, "strong.pickle")
    sequence_path = os.path.join(data_dir, "experiments.json")

    unit_cell = uctbx.unit_cell((58.373, 58.373, 155.939, 90, 90, 120))
    hall_symbol = '-R 3 2"'

    extra_args = [
        "bin_size_fraction=0.25",
        'known_symmetry.space_group="Hall: %s"' % hall_symbol.replace('"', '\\"'),
    ]
    if specify_unit_cell:
        extra_args.append(
            'known_symmetry.unit_cell="%s %s %s %s %s %s"' % unit_cell.parameters()
        )

    expected_unit_cell = unit_cell
    expected_hall_symbol = hall_symbol
    expected_rmsds = (0.08, 0.11, 0.004)

    run_indexing(
        pickle_path,
        sequence_path,
        tmpdir,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )


@pytest.fixture(scope="session")
def insulin_spotfinding(dials_data, tmpdir_factory):
    """Return experiment and reflection files for 2 images of the insulin dataset"""

    data_dir = dials_data("insulin")
    tmpdir = tmpdir_factory.mktemp("insulin")

    command = ["dials.import"]
    for i, image_path in enumerate(("insulin_1_001.img", "insulin_1_045.img")):
        target = "image_00%i.img" % (i + 1)
        data_dir.join(image_path).copy(tmpdir.join(target))
        command.append(target)

    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr

    experiment = tmpdir.join("imported.expt")
    assert experiment.check()

    command = ["dials.find_spots", experiment]
    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr

    reflections = tmpdir.join("strong.refl")
    assert reflections.check()

    return experiment, reflections


@pytest.mark.parametrize("method", ["fft3d", "fft1d", "real_space_grid_search"])
def test_index_insulin_multi_sequence(insulin_spotfinding, tmpdir, method):
    experiment, reflections = insulin_spotfinding
    expected_unit_cell = uctbx.unit_cell(
        (78.163, 78.163, 78.163, 90.000, 90.000, 90.000)
    )
    expected_hall_symbol = " I 2 2 3"
    expected_rmsds = (0.05, 0.06, 0.01)
    extra_args = [
        'known_symmetry.unit_cell="%s %s %s %s %s %s"'
        % expected_unit_cell.parameters(),
        'known_symmetry.space_group="Hall: %s"' % expected_hall_symbol,
        "indexing.method=%s" % method,
        "treat_single_image_as_still=False",
    ]
    run_indexing(
        reflections,
        experiment,
        tmpdir,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )


@pytest.fixture(scope="session")
def insulin_spotfinding_stills(dials_data, tmpdir_factory):
    """Return experiment and reflection files for 1 image of the insulin
    dataset treated as still image"""

    data_dir = dials_data("insulin")
    tmpdir = tmpdir_factory.mktemp("insulin")

    command = [
        "dials.import",
        "convert_sequences_to_stills=True",
        data_dir.join("insulin_1_001.img"),
    ]
    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr

    experiment = tmpdir.join("imported.expt")
    assert experiment.check()

    command = ["dials.find_spots", experiment]
    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr

    reflections = tmpdir.join("strong.refl")
    assert reflections.check()

    return experiment, reflections


@pytest.mark.parametrize("method", ["fft3d", "fft1d", "real_space_grid_search"])
def test_index_insulin_force_stills(insulin_spotfinding_stills, tmpdir, method):
    experiment, reflections = insulin_spotfinding_stills
    expected_unit_cell = uctbx.unit_cell(
        (78.163, 78.163, 78.163, 90.000, 90.000, 90.000)
    )
    expected_hall_symbol = " I 2 2 3"
    expected_rmsds = (0.05, 0.06, 0.01)

    extra_args = [
        "stills.indexer=stills",
        'known_symmetry.unit_cell="%s %s %s %s %s %s"'
        % expected_unit_cell.parameters(),
        'known_symmetry.space_group="Hall: %s"' % expected_hall_symbol,
        "indexing.method=%s" % method,
    ]

    run_indexing(
        reflections,
        experiment,
        tmpdir,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )


def test_multiple_experiments(dials_regression, tmpdir):
    # Test indexing 4 lysozyme still shots in a single dials.index job
    #   - the first image doesn't index
    #   - the last three images do index
    data_dir = os.path.join(
        dials_regression, "indexing_test_data", "i24_lysozyme_stills"
    )
    pickle_path = os.path.join(data_dir, "strong.pickle")
    experiments_json = os.path.join(data_dir, "imported_experiments.json")

    expected_unit_cell = uctbx.unit_cell((38.06, 78.78, 78.91, 90, 90, 90))
    expected_hall_symbol = " P 1"
    expected_rmsds = (0.1, 0.07, 0.0)

    extra_args = [
        "stills.indexer=sequences",
        "joint_indexing=False",
        "outlier.algorithm=sauter_poon",
    ]

    run_indexing(
        pickle_path,
        experiments_json,
        tmpdir,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
        n_expected_lattices=3,
        relative_length_tolerance=0.01,
    )


def test_index_4rotation(dials_regression, tmpdir):
    data_dir = os.path.join(dials_regression, "indexing_test_data", "4rotation")
    pickle_path = os.path.join(data_dir, "strong.pickle")
    sequence_path = os.path.join(data_dir, "experiments.json")
    extra_args = [
        "max_refine=10",
        "reflections_per_degree=50",
        "known_symmetry.space_group=R3",
        "n_macro_cycles=3",
    ]
    expected_unit_cell = uctbx.unit_cell((48.397, 48.397, 284.767, 90, 90, 120))
    expected_rmsds = (0.06, 0.08, 0.22)
    expected_hall_symbol = " R 3"

    result = run_indexing(
        pickle_path,
        sequence_path,
        tmpdir,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )
    assert len(result.indexed_reflections) > 276800, len(result.indexed_reflections)


def test_index_small_molecule_multi_sequence_4(dials_regression, tmpdir):
    # test for small molecule multi-sequence indexing, 4 sequences with different values
    # of goniometer.fixed_rotation()
    data_dir = os.path.join(dials_regression, "indexing_test_data", "multi_sweep")
    pickle_paths = [
        glob.glob(
            os.path.join(data_dir, "SWEEP%i" % (i + 1), "index", "*_strong.pickle")
        )[0]
        for i in range(4)
    ]
    sequence_paths = [
        glob.glob(
            os.path.join(data_dir, "SWEEP%i" % (i + 1), "index", "experiments.json")
        )[0]
        for i in range(4)
    ]
    extra_args = ["known_symmetry.space_group=I4", "filter_ice=False"]
    expected_unit_cell = uctbx.unit_cell((7.310, 7.310, 6.820, 90.000, 90.000, 90.000))
    expected_rmsds = (0.10, 0.7, 0.5)
    expected_hall_symbol = " I 4"

    result = run_indexing(
        pickle_paths,
        sequence_paths,
        tmpdir,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )
    assert len(result.indexed_reflections) > 1250, len(result.indexed_reflections)


def test_index_small_molecule_multi_sequence_3(dials_regression, tmpdir):
    # test for small molecule multi-sequence indexing, 3 sequences with different values
    # of goniometer setting rotation (i.e. phi scans)
    data_dir = os.path.join(dials_regression, "dials-191")
    pickle_paths = [
        glob.glob(os.path.join(data_dir, "*SWEEP%i*_strong.pickle" % (i + 1)))[0]
        for i in range(3)
    ]
    sequence_paths = [
        glob.glob(os.path.join(data_dir, "*SWEEP%i*_experiments.json" % (i + 1)))[0]
        for i in range(3)
    ]
    extra_args = ["filter_ice=False"]
    expected_unit_cell = uctbx.unit_cell(
        (9.440, 15.313, 17.126, 90.073, 90.106, 79.248)
    )
    expected_rmsds = (0.32, 0.34, 0.005)
    expected_hall_symbol = " P 1"

    result = run_indexing(
        pickle_paths,
        sequence_paths,
        tmpdir,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )
    assert len(result.indexed_reflections) > 12000, len(result.indexed_reflections)
    # expect at least indexed 2000 reflections per experiment
    for i in range(3):
        assert (result.indexed_reflections["id"] == i).count(True) > 2000


def test_index_small_molecule_ice_max_cell(dials_regression, tmpdir):
    # test for small molecule indexing: presence of ice rings makes max-cell
    # estimation tricky
    data_dir = os.path.join(dials_regression, "indexing_test_data", "MXSW-904")
    pickle_path = os.path.join(data_dir, "1_SWEEP1_strong.pickle")
    experiments = os.path.join(data_dir, "1_SWEEP1_experiments.json")
    extra_args = ["filter_ice=False"]
    expected_unit_cell = uctbx.unit_cell((11.72, 11.72, 11.74, 109.08, 109.24, 108.99))
    expected_rmsds = (0.06, 0.05, 0.04)
    expected_hall_symbol = " P 1"

    result = run_indexing(
        pickle_path,
        experiments,
        tmpdir,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )
    assert len(result.indexed_reflections) > 1300, len(result.indexed_reflections)


def test_refinement_failure_on_max_lattices_a15(dials_regression, tmpdir):
    """Problem: Sometimes there is enough data to index, but not enough to
    refine. If this happens in the (N>1)th crystal of max_lattices, then
    all existing solutions are also dropped."""
    data_dir = os.path.join(dials_regression, "indexing_test_data", "lattice_failures")

    result = procrunner.run(
        [
            "dials.index",
            os.path.join(data_dir, "lpe4-2-a15_strong.pickle"),
            os.path.join(data_dir, "lpe4-2-a15_datablock.json"),
            "max_lattices=3",
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("indexed.refl").check() and tmpdir.join("indexed.expt").check()
    experiments_list = load.experiment_list(
        tmpdir.join("indexed.expt").strpath, check_format=False
    )
    assert len(experiments_list) == 2

    # now try to reindex with existing model
    result = procrunner.run(
        [
            "dials.index",
            tmpdir.join("indexed.expt").strpath,
            os.path.join(data_dir, "lpe4-2-a15_strong.pickle"),
            "max_lattices=2",
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("indexed.refl").check() and tmpdir.join("indexed.expt").check()
    experiments_list = load.experiment_list(
        tmpdir.join("indexed.expt").strpath, check_format=False
    )
    assert len(experiments_list) == 2


def test_stills_indexer_multi_lattice_bug_MosaicSauter2014(dials_regression, tmpdir):
    """ Problem: In stills_indexer, before calling the refine function, the experiment list contains a list of
        dxtbx crystal models (that are not MosaicSauter2014 models). The conversion to MosaicSauter2014 is made
        during the refine step when functions from nave_parameters is called. If the experiment list contains
        more than 1 experiment, for eg. multiple lattices, only the first crystal gets assigned mosaicity. In
        actuality, all crystal models should be assigned mosaicity. This test only compares whether or not all crystal models
        have been assigned a MosaicSauter2014 model.  """

    from dxtbx.model.experiment_list import ExperimentListFactory
    from dxtbx.model.experiment_list import Experiment, ExperimentList
    from dials.array_family import flex
    from dxtbx.model import Crystal
    from dials.algorithms.indexing.stills_indexer import StillsIndexer
    from dials.command_line.stills_process import (
        phil_scope as stills_process_phil_scope,
    )
    import dxtbx_model_ext  # needed for comparison of types

    experiment_data = os.path.join(
        dials_regression,
        "refinement_test_data",
        "cspad_refinement",
        "cspad_refined_experiments_step6_level2_300.json",
    )
    reflection_data = os.path.join(
        dials_regression,
        "refinement_test_data",
        "cspad_refinement",
        "cspad_reflections_step7_300.pickle",
    )

    refl = flex.reflection_table.from_file(reflection_data)
    explist = ExperimentListFactory.from_json_file(experiment_data, check_format=False)[
        0:2
    ]
    reflist = refl.select(refl["id"] < 2)  # Only use the first 2 for convenience
    # Construct crystal models that don't have mosaicity. These A,B,C values are the same
    # as read in from the dials_regression folder
    # Crystal-0
    cs0 = Crystal(explist[0].crystal)
    exp0 = Experiment(
        imageset=explist[0].imageset,
        beam=explist[0].beam,
        detector=explist[0].detector,
        goniometer=None,
        scan=None,
        crystal=cs0,
    )

    # Crystal-1
    cs1 = Crystal(explist[1].crystal)
    exp1 = Experiment(
        imageset=explist[1].imageset,
        beam=explist[1].beam,
        detector=explist[1].detector,
        goniometer=None,
        scan=None,
        crystal=cs1,
    )
    # Construct a new experiment_list that will be passed on for refinement
    unrefined_explist = ExperimentList([exp0, exp1])
    # Get default params from stills_process and construct StillsIndexer, then run refinement
    params = stills_process_phil_scope.extract()
    SI = StillsIndexer(reflist, unrefined_explist, params=params)
    refined_explist, new_reflist = SI.refine(unrefined_explist, reflist)
    # Now check whether the models have mosaicity after stills_indexer refinement
    # Also check that mosaicity values are within expected limits
    for ii, crys in enumerate(refined_explist.crystals()):
        assert isinstance(crys, dxtbx_model_ext.MosaicCrystalSauter2014)
        if ii == 0:
            assert crys.get_domain_size_ang() == pytest.approx(2242.0, rel=0.1)
        if ii == 1:
            assert crys.get_domain_size_ang() == pytest.approx(2689.0, rel=0.1)


@pytest.mark.parametrize(
    "indexer_type,fix_cell", (("sequences", False), ("stills", True))
)
def test_index_ED_still_low_res_spot_match(dials_data, tmpdir, indexer_type, fix_cell):

    # test indexing from a single simulated lysozyme ED still

    image_path = dials_data("image_examples").join("simtbx_FormatSMVJHSim_001.img")

    command = ["dials.import", image_path]
    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr

    experiment = tmpdir.join("imported.expt")
    assert experiment.check()

    command = ["dials.find_spots", experiment]
    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr

    reflections = tmpdir.join("strong.refl")

    extra_args = [
        "indexing.method=low_res_spot_match",
        "known_symmetry.space_group=P43212",
        "known_symmetry.unit_cell=78.84,78.84,38.29,90,90,90",
        "stills.indexer=" + indexer_type,
        "n_macro_cycles=2",
        "detector.fix_list=Dist",
    ]
    if fix_cell:
        extra_args += ["crystal.fix=cell"]

    expected_unit_cell = uctbx.unit_cell((78.84, 78.84, 38.29, 90, 90, 90))
    expected_rmsds = (0.007, 0.006, 0.000)
    expected_hall_symbol = " P 4nw 2abw"

    run_indexing(
        reflections,
        experiment,
        tmpdir,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )


def test_real_space_grid_search_no_unit_cell(dials_regression, tmpdir):
    data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")
    experiments_json = os.path.join(data_dir, "experiments_import.json")
    pickle_path = os.path.join(data_dir, "full.pickle")
    commands = [
        "dials.index",
        experiments_json,
        pickle_path,
        "indexing.method=real_space_grid_search",
    ]
    result = procrunner.run(commands, working_directory=tmpdir)
    assert result.stderr
    assert (
        result.stderr.strip()
        == b"Target unit cell must be provided for real_space_grid_search"
    )


def test_index_known_orientation(dials_data, tmpdir):
    data_dir = dials_data("vmxi_proteinase_k_sweeps")
    experiments_json = data_dir.join("experiments_0.json").strpath
    reflections = data_dir.join("reflections_0.pickle").strpath

    expected_unit_cell = uctbx.unit_cell((68.395, 68.395, 104, 90, 90, 90))
    expected_rmsds = (0.013, 0.012, 0.008)
    expected_hall_symbol = " P 4"

    run_indexing(
        reflections,
        experiments_json,
        tmpdir,
        [],
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )


def test_all_expt_ids_have_expts(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.index",
            dials_data("vmxi_thaumatin_grid_index").join("split_07602.expt"),
            dials_data("vmxi_thaumatin_grid_index").join("split_07602.refl"),
            "stills.indexer=sequences",
            "indexing.method=real_space_grid_search",
            "space_group=P4",
            "unit_cell=58,58,150,90,90,90",
            "max_lattices=8",
            "beam.fix=all",
            "detector.fix=all",
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("indexed.expt").check(file=1)
    assert tmpdir.join("indexed.refl").check(file=1)

    refl = flex.reflection_table.from_file(tmpdir / "indexed.refl")
    expt = ExperimentList.from_file(tmpdir / "indexed.expt")

    assert flex.max(refl["id"]) + 1 == len(expt)
