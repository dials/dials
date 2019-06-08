from __future__ import absolute_import, division, print_function

import glob
import os
import pytest

from libtbx import easy_run
from scitbx import matrix
from cctbx import uctbx
from dxtbx.serialize import load
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


class RunOneIndexing(object):
    def __init__(
        self,
        pickle_path,
        sweep_path,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
        n_expected_lattices=1,
        relative_length_tolerance=0.005,
        absolute_angle_tolerance=0.5,
    ):

        args = ["dials.index", pickle_path, sweep_path] + extra_args
        command = " ".join(args)
        print(command)
        easy_run.fully_buffered(command=command).raise_if_errors()
        assert os.path.exists("indexed_experiments.json")
        experiments_list = load.experiment_list(
            "indexed_experiments.json", check_format=False
        )
        assert len(experiments_list.crystals()) == n_expected_lattices, (
            len(experiments_list.crystals()),
            n_expected_lattices,
        )
        assert os.path.exists("indexed.pickle")

        self.indexed_reflections = flex.reflection_table.from_file("indexed.pickle")

        self.experiments = experiments_list
        for i in range(len(experiments_list)):
            experiment = experiments_list[i]
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
            reflections = self.indexed_reflections.select(
                self.indexed_reflections["id"] == i
            )
            mi = reflections["miller_index"]
            assert (mi != (0, 0, 0)).count(False) == 0
            reflections = reflections.select(mi != (0, 0, 0))
            self.rmsds = self.get_rmsds_obs_pred(reflections, experiment)
            for actual, expected in zip(self.rmsds, expected_rmsds):
                assert actual <= expected, "%s %s" % (self.rmsds, expected_rmsds)

    def get_rmsds_obs_pred(self, observations, experiment):
        reflections = observations.select(
            observations.get_flags(observations.flags.used_in_refinement)
        )
        assert len(reflections) > 0
        obs_x, obs_y, obs_z = reflections["xyzobs.mm.value"].parts()
        calc_x, calc_y, calc_z = reflections["xyzcal.mm"].parts()
        rmsd_x = flex.mean(flex.pow2(obs_x - calc_x)) ** 0.5
        rmsd_y = flex.mean(flex.pow2(obs_y - calc_y)) ** 0.5
        rmsd_z = flex.mean(flex.pow2(obs_z - calc_z)) ** 0.5
        return (rmsd_x, rmsd_y, rmsd_z)


def test_index_i04_weak_data_fft3d(dials_regression, tmpdir):
    # thaumatin
    data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")
    pickle_path = os.path.join(data_dir, "full.pickle")
    sweep_path = os.path.join(data_dir, "experiments_import.json")
    extra_args = [
        "bin_size_fraction=0.25",
        "image_range=1,20",
        "image_range=250,270",
        "image_range=520,540",
    ]
    expected_unit_cell = uctbx.unit_cell((57.7, 57.7, 149.8, 90, 90, 90))
    expected_rmsds = (0.05, 0.04, 0.0005)
    expected_hall_symbol = " P 1"

    with tmpdir.as_cwd():
        RunOneIndexing(
            pickle_path,
            sweep_path,
            extra_args,
            expected_unit_cell,
            expected_rmsds,
            expected_hall_symbol,
        )


@pytest.mark.skip(reason="cluster_analysis_search=True not implemented")
def test_index_cluster_analysis_search(dials_regression, tmpdir):
    pytest.importorskip("scipy")
    pytest.importorskip("sklearn")
    pytest.importorskip("networkx")

    # thaumatin
    data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")
    pickle_path = os.path.join(data_dir, "full.pickle")
    sweep_path = os.path.join(data_dir, "experiments_import.json")
    extra_args = [
        "cluster_analysis_search=True",
        "n_macro_cycles=3",
        "bin_size_fraction=0.25",
        "reciprocal_space_grid.d_min=4",
    ]
    expected_unit_cell = uctbx.unit_cell((58, 58, 150, 90, 90, 90))
    expected_rmsds = (0.05, 0.04, 0.0004)
    expected_hall_symbol = " P 1"

    with tmpdir.as_cwd():
        RunOneIndexing(
            pickle_path,
            sweep_path,
            extra_args,
            expected_unit_cell,
            expected_rmsds,
            expected_hall_symbol,
        )


@pytest.mark.skip(reason="cluster_analysis_search=True not implemented")
def test_index_cluster_analysis_search_with_symmetry(dials_regression, tmpdir):
    pytest.importorskip("scipy")
    pytest.importorskip("sklearn")
    pytest.importorskip("networkx")

    # thaumatin
    data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")
    pickle_path = os.path.join(data_dir, "full.pickle")
    sweep_path = os.path.join(data_dir, "experiments_import.json")
    extra_args = [
        "cluster_analysis_search=True",
        "n_macro_cycles=3",
        "bin_size_fraction=0.25",
        "reciprocal_space_grid.d_min=4",
    ]
    expected_unit_cell = uctbx.unit_cell((58, 58, 150, 90, 90, 90))
    expected_rmsds = (0.05, 0.042, 0.0004)

    # now enforce symmetry
    extra_args.append("known_symmetry.space_group=P4")
    expected_hall_symbol = " P 4"

    with tmpdir.as_cwd():
        result = RunOneIndexing(
            pickle_path,
            sweep_path,
            extra_args,
            expected_unit_cell,
            expected_rmsds,
            expected_hall_symbol,
        )

    a, b, c = map(matrix.col, result.experiments[0].crystal.get_real_space_vectors())
    assert a.length() == pytest.approx(b.length())
    assert c.length() > b.length()
    assert a.angle(b, deg=True) == pytest.approx(90)
    assert b.angle(c, deg=True) == pytest.approx(90)
    assert c.angle(a, deg=True) == pytest.approx(90)


@pytest.mark.skip(reason="cluster_analysis_search=True not implemented")
def test_index_trypsin_single_lattice(dials_regression, tmpdir):
    pytest.importorskip("scipy")
    pytest.importorskip("sklearn")
    pytest.importorskip("networkx")

    # trypsin
    data_dir = os.path.join(dials_regression, "indexing_test_data", "trypsin")
    pickle_path = os.path.join(data_dir, "P1_X6_1.pickle")
    sweep_path = os.path.join(data_dir, "experiments_P1_X6_1.json")
    extra_args = [
        "cluster_analysis_search=True",
        "n_macro_cycles=3",
        "reciprocal_space_grid.d_min=4",
        "filter_overlaps=False",  # P1_X6_1.pickle does not contain bbox!
        "image_range=0,50",
        "image_range=450,500",
        "image_range=850,900",
    ]
    expected_unit_cell = uctbx.unit_cell((54.3, 58.3, 66.5, 90, 90, 90))
    expected_rmsds = (0.061, 0.06, 0.00042)
    expected_hall_symbol = " P 1"

    with tmpdir.as_cwd():
        RunOneIndexing(
            pickle_path,
            sweep_path,
            extra_args,
            expected_unit_cell,
            expected_rmsds,
            expected_hall_symbol,
        )


@pytest.mark.skip(reason="cluster_analysis_search=True not implemented")
def test_index_trypsin_two_lattice(dials_regression, tmpdir):
    pytest.importorskip("scipy")
    pytest.importorskip("sklearn")
    pytest.importorskip("networkx")

    # synthetic trypsin multi-lattice dataset (2 lattices)
    data_dir = os.path.join(dials_regression, "indexing_test_data", "trypsin")
    pickle_path = os.path.join(data_dir, "P1_X6_1_2.pickle")
    sweep_path = os.path.join(data_dir, "experiments_P1_X6_1_2.json")
    extra_args = [
        "cluster_analysis_search=True",
        "reflections_per_degree=10",
        "n_macro_cycles=2",
        "reciprocal_space_grid.d_min=4",
        "max_cell=70",
        "image_range=0,50",
        "image_range=450,500",
        "image_range=850,900",
        "max_lattices=2",
    ]
    expected_unit_cell = uctbx.unit_cell((54.3, 58.3, 66.5, 90, 90, 90))
    expected_rmsds = (0.1, 0.1, 0.005)
    expected_hall_symbol = " P 1"
    n_expected_lattices = 2

    with tmpdir.as_cwd():
        RunOneIndexing(
            pickle_path,
            sweep_path,
            extra_args,
            expected_unit_cell,
            expected_rmsds,
            expected_hall_symbol,
            n_expected_lattices=n_expected_lattices,
            relative_length_tolerance=0.02,
            absolute_angle_tolerance=1,
        )


@pytest.mark.skip(reason="cluster_analysis_search=True not implemented")
def test_index_trypsin_three_lattice_cluster_analysis_search(dials_regression, tmpdir):
    pytest.importorskip("scipy")
    pytest.importorskip("sklearn")
    pytest.importorskip("networkx")

    # synthetic trypsin multi-lattice dataset (3 lattices)
    data_dir = os.path.join(dials_regression, "indexing_test_data", "trypsin")
    pickle_path = os.path.join(data_dir, "P1_X6_1_2_3.pickle")
    sweep_path = os.path.join(data_dir, "experiments_P1_X6_1_2_3.json")
    extra_args = [
        "cluster_analysis_search=True",
        "reflections_per_degree=10",
        "n_macro_cycles=2",
        "reciprocal_space_grid.d_min=4",
        "max_cell=70",
        "max_lattices=3",
    ]
    expected_unit_cell = uctbx.unit_cell((54.3, 58.3, 66.5, 90, 90, 90))
    expected_rmsds = (0.24, 0.30, 0.006)
    expected_hall_symbol = " P 1"
    n_expected_lattices = 3

    with tmpdir.as_cwd():
        RunOneIndexing(
            pickle_path,
            sweep_path,
            extra_args,
            expected_unit_cell,
            expected_rmsds,
            expected_hall_symbol,
            n_expected_lattices=n_expected_lattices,
            relative_length_tolerance=0.01,
            absolute_angle_tolerance=1,
        )


@pytest.mark.skip(reason="cluster_analysis_search=True not implemented")
def test_index_trypsin_four_lattice_P1(dials_regression, tmpdir):
    pytest.importorskip("scipy")
    pytest.importorskip("sklearn")
    pytest.importorskip("networkx")

    # synthetic trypsin multi-lattice dataset (4 lattices)
    data_dir = os.path.join(dials_regression, "indexing_test_data", "trypsin")
    pickle_path = os.path.join(data_dir, "P1_X6_1_2_3_4.pickle")
    sweep_path = os.path.join(data_dir, "experiments_P1_X6_1_2_3_4.json")
    extra_args = [
        "cluster_analysis_search=True",
        "reflections_per_degree=10",
        "n_macro_cycles=2",
        "reciprocal_space_grid.d_min=4",
        "max_cell=70",
        "max_lattices=4",
    ]
    expected_unit_cell = uctbx.unit_cell((54.3, 58.3, 66.5, 90, 90, 90))
    expected_rmsds = (0.24, 0.23, 0.006)
    expected_hall_symbol = " P 1"
    n_expected_lattices = 4

    with tmpdir.as_cwd():
        RunOneIndexing(
            pickle_path,
            sweep_path,
            extra_args,
            expected_unit_cell,
            expected_rmsds,
            expected_hall_symbol,
            n_expected_lattices=n_expected_lattices,
            relative_length_tolerance=0.01,
            absolute_angle_tolerance=1,
        )


def test_index_trypsin_four_lattice_P212121(dials_regression, tmpdir):
    # synthetic trypsin multi-lattice dataset (4 lattices)
    data_dir = os.path.join(dials_regression, "indexing_test_data", "trypsin")
    pickle_path = os.path.join(data_dir, "P1_X6_1_2_3_4.pickle")
    sweep_path = os.path.join(data_dir, "experiments_P1_X6_1_2_3_4.json")
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

    with tmpdir.as_cwd():
        RunOneIndexing(
            pickle_path,
            sweep_path,
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
    sweep_path = os.path.join(data_dir, "experiments_import.json")
    extra_args = [
        "n_macro_cycles=2",
        "indexing.method=fft1d",
        "bin_size_fraction=0.25",
        "image_range=1,20",
        "image_range=250,270",
        "image_range=520,540",
    ]
    expected_unit_cell = uctbx.unit_cell((58, 58, 150, 90, 90, 90))
    expected_rmsds = (0.06, 0.05, 0.0005)
    expected_hall_symbol = " P 1"

    with tmpdir.as_cwd():
        RunOneIndexing(
            pickle_path,
            sweep_path,
            extra_args,
            expected_unit_cell,
            expected_rmsds,
            expected_hall_symbol,
        )


def test_index_trypsin_index_assignment_local(dials_regression, tmpdir):
    # synthetic trypsin multi-lattice dataset (3 lattices)
    data_dir = os.path.join(dials_regression, "indexing_test_data", "trypsin")
    pickle_path = os.path.join(data_dir, "P1_X6_1_2_3.pickle")
    sweep_path = os.path.join(data_dir, "experiments_P1_X6_1_2_3.json")
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

    with tmpdir.as_cwd():
        RunOneIndexing(
            pickle_path,
            sweep_path,
            extra_args,
            expected_unit_cell,
            expected_rmsds,
            expected_hall_symbol,
            n_expected_lattices=n_expected_lattices,
            relative_length_tolerance=0.02,
            absolute_angle_tolerance=1,
        )


def test_index_peak_search_clean(dials_regression, tmpdir):
    pytest.importorskip("scipy")
    pytest.importorskip("sklearn")
    pytest.importorskip("networkx")

    # test indexing from single image of i04_weak_data
    data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")
    pickle_path = os.path.join(data_dir, "first_image.pickle")
    sweep_path = os.path.join(data_dir, "experiments_import.json")
    extra_args = [
        "indexing.method=fft3d",
        "known_symmetry.space_group=P4",
        "known_symmetry.unit_cell=57.8,57.8,150,90,90,90",
        "peak_search=clean",
        "cluster_analysis_search=True",
        "min_samples=15",
        "n_macro_cycles=4",
        "reciprocal_space_grid.d_min=4",
    ]

    expected_unit_cell = uctbx.unit_cell((57.8, 57.8, 150, 90, 90, 90))
    expected_rmsds = (0.06, 0.07, 0.003)
    expected_hall_symbol = " P 4"

    with tmpdir.as_cwd():
        RunOneIndexing(
            pickle_path,
            sweep_path,
            extra_args,
            expected_unit_cell,
            expected_rmsds,
            expected_hall_symbol,
        )


def test_index_imosflm_tutorial(dials_regression, tmpdir):
    # test on spots derived from imosflm tutorial data:
    # http://www.ccp4.ac.uk/courses/BCA2005/tutorials/dataproc-tutorial.html
    data_dir = os.path.join(dials_regression, "indexing_test_data", "imosflm_hg_mar")
    pickle_path = os.path.join(data_dir, "strong.pickle")
    sweep_path = os.path.join(data_dir, "experiments.json")

    unit_cell = uctbx.unit_cell((58.373, 58.373, 155.939, 90, 90, 120))
    hall_symbol = '-R 3 2"'

    for uc, hall in ((unit_cell, hall_symbol), (None, hall_symbol)):
        extra_args = ["bin_size_fraction=0.25"]
        if uc is not None:
            extra_args.append(
                'known_symmetry.unit_cell="%s %s %s %s %s %s"' % unit_cell.parameters()
            )
        if hall is not None:
            extra_args.append(
                'known_symmetry.space_group="Hall: %s"' % hall.replace('"', '\\"')
            )

        expected_unit_cell = unit_cell
        if hall is not None:
            expected_hall_symbol = hall
        else:
            expected_hall_symbol = " P 1"
        expected_rmsds = (0.08, 0.11, 0.004)

        with tmpdir.as_cwd():
            RunOneIndexing(
                pickle_path,
                sweep_path,
                extra_args,
                expected_unit_cell,
                expected_rmsds,
                expected_hall_symbol,
            )


@pytest.mark.parametrize("method", ["fft3d", "fft1d", "real_space_grid_search"])
def test_index_insulin_multi_sweep(dials_data, run_in_tmpdir, method):
    data_dir = dials_data("insulin")

    args = ["dials.import", "allow_multiple_sweeps=True"]
    for i, image_path in enumerate(("insulin_1_001.img", "insulin_1_045.img")):
        target = "image_00%i.img" % (i + 1)
        data_dir.join(image_path).copy(run_in_tmpdir.join(target))
        args.append(target)

    command = " ".join(args)
    # print(command)
    easy_run.fully_buffered(command=command).raise_if_errors()

    experiments_json = "imported_experiments.json"

    args = ["dials.find_spots", experiments_json]

    command = " ".join(args)
    print(command)
    easy_run.fully_buffered(command=command).raise_if_errors()
    pickle_path = "strong.pickle"
    assert os.path.exists(pickle_path)

    expected_unit_cell = uctbx.unit_cell(
        (78.163, 78.163, 78.163, 90.000, 90.000, 90.000)
    )
    expected_hall_symbol = " I 2 2 3"
    expected_rmsds = (0.05, 0.06, 0.01)

    extra_args = []
    extra_args.append(
        'known_symmetry.unit_cell="%s %s %s %s %s %s"' % expected_unit_cell.parameters()
    )
    extra_args.append('known_symmetry.space_group="Hall: %s"' % expected_hall_symbol)
    extra_args.append("indexing.method=%s" % method)
    extra_args.append("treat_single_image_as_still=False")

    RunOneIndexing(
        pickle_path,
        experiments_json,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )


@pytest.mark.parametrize("method", ["fft3d", "fft1d", "real_space_grid_search"])
def test_index_insulin_force_stills(dials_data, run_in_tmpdir, method):
    data_dir = dials_data("insulin")

    args = [
        "dials.import",
        "convert_sweeps_to_stills=True",
        data_dir.join("insulin_1_001.img").strpath,
    ]
    command = " ".join(args)
    # print(command)
    easy_run.fully_buffered(command=command).raise_if_errors()

    experiments_json = "imported_experiments.json"
    args = ["dials.find_spots", experiments_json]
    command = " ".join(args)
    easy_run.fully_buffered(command=command).raise_if_errors()
    pickle_path = "strong.pickle"
    assert os.path.exists(pickle_path)

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

    RunOneIndexing(
        pickle_path,
        experiments_json,
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

    extra_args = ["stills.indexer=sweeps", "joint_indexing=False"]

    with tmpdir.as_cwd():
        RunOneIndexing(
            pickle_path,
            experiments_json,
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
    sweep_path = os.path.join(data_dir, "experiments.json")
    extra_args = [
        "max_refine=10",
        "reflections_per_degree=50",
        "known_symmetry.space_group=R3",
        "n_macro_cycles=3",
    ]
    expected_unit_cell = uctbx.unit_cell((48.397, 48.397, 284.767, 90, 90, 120))
    expected_rmsds = (0.06, 0.08, 0.22)
    expected_hall_symbol = " R 3"

    with tmpdir.as_cwd():
        result = RunOneIndexing(
            pickle_path,
            sweep_path,
            extra_args,
            expected_unit_cell,
            expected_rmsds,
            expected_hall_symbol,
        )
        assert len(result.indexed_reflections) > 276800, len(result.indexed_reflections)


def test_index_small_molecule_multi_sweep_4(dials_regression, tmpdir):
    # test for small molecule multi-sweep indexing, 4 sweeps with different values
    # of goniometer.fixed_rotation()
    data_dir = os.path.join(dials_regression, "indexing_test_data", "multi_sweep")
    pickle_paths = [
        glob.glob(
            os.path.join(data_dir, "SWEEP%i" % (i + 1), "index", "*_strong.pickle")
        )[0]
        for i in range(4)
    ]
    sweep_paths = [
        glob.glob(
            os.path.join(data_dir, "SWEEP%i" % (i + 1), "index", "experiments.json")
        )[0]
        for i in range(4)
    ]
    extra_args = ["known_symmetry.space_group=I4", "filter_ice=False"]
    expected_unit_cell = uctbx.unit_cell((7.310, 7.310, 6.820, 90.000, 90.000, 90.000))
    expected_rmsds = (0.10, 0.7, 0.5)
    expected_hall_symbol = " I 4"

    with tmpdir.as_cwd():
        result = RunOneIndexing(
            " ".join(pickle_paths),
            " ".join(sweep_paths),
            extra_args,
            expected_unit_cell,
            expected_rmsds,
            expected_hall_symbol,
        )
        assert len(result.indexed_reflections) > 1250, len(result.indexed_reflections)


def test_index_small_molecule_multi_sweep_3(dials_regression, tmpdir):
    # test for small molecule multi-sweep indexing, 3 sweeps with different values
    # of goniometer setting rotation (i.e. phi scans)
    data_dir = os.path.join(dials_regression, "dials-191")
    pickle_paths = [
        glob.glob(os.path.join(data_dir, "*SWEEP%i*_strong.pickle" % (i + 1)))[0]
        for i in range(3)
    ]
    sweep_paths = [
        glob.glob(os.path.join(data_dir, "*SWEEP%i*_experiments.json" % (i + 1)))[0]
        for i in range(3)
    ]
    extra_args = ["filter_ice=False"]
    expected_unit_cell = uctbx.unit_cell(
        (9.440, 15.313, 17.126, 90.073, 90.106, 79.248)
    )
    expected_rmsds = (0.32, 0.34, 0.005)
    expected_hall_symbol = " P 1"

    with tmpdir.as_cwd():
        result = RunOneIndexing(
            " ".join(pickle_paths),
            " ".join(sweep_paths),
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

    with tmpdir.as_cwd():
        result = RunOneIndexing(
            pickle_path,
            experiments,
            extra_args,
            expected_unit_cell,
            expected_rmsds,
            expected_hall_symbol,
        )
        assert len(result.indexed_reflections) > 1300, len(result.indexed_reflections)


def test_refinement_failure_on_max_lattices_a15(dials_regression, run_in_tmpdir):
    """Problem: Sometimes there is enough data to index, but not enough to
    refine. If this happens in the (N>1)th crystal of max_lattices, then
    all existing solutions are also dropped."""
    data_dir = os.path.join(dials_regression, "indexing_test_data", "lattice_failures")

    cmd = [
        "dials.index",
        os.path.join(data_dir, "lpe4-2-a15_strong.pickle"),
        os.path.join(data_dir, "lpe4-2-a15_datablock.json"),
        "max_lattices=3",
    ]
    easy_run.fully_buffered(command=" ".join(cmd)).raise_if_errors()
    assert os.path.isfile("indexed.pickle") and os.path.isfile(
        "indexed_experiments.json"
    )
    experiments_list = load.experiment_list(
        "indexed_experiments.json", check_format=False
    )
    assert len(experiments_list) == 2
