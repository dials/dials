"""
Test combination of multiple experiments and reflections files.
"""

from __future__ import annotations

import copy
import os
import shutil
import subprocess

import pytest

from dxtbx.model.experiment_list import ExperimentListFactory
from dxtbx.serialize import load

from dials.array_family import flex
from dials.command_line.combine_experiments import (
    combine_experiments,
    combine_experiments_no_reflections,
    phil_scope,
)


def test(dials_data, tmp_path):
    data_dir = dials_data("polyhedra_narrow_wedges", pathlib=True)

    input_range = list(range(2, 49))
    for i in (8, 10, 15, 16, 34, 39, 45):
        input_range.remove(i)

    phil_input = (
        "\n".join(
            (
                f"  input.experiments={data_dir}/sweep_{i:03d}_experiments.json\n"
                + f"  input.reflections={data_dir}/sweep_{i:03d}_reflections.pickle"
            )
            for i in input_range
        )
        + """
 reference_from_experiment.beam=0
 reference_from_experiment.scan=0
 reference_from_experiment.goniometer=0
 reference_from_experiment.detector=0
 """
    )
    tmp_path.joinpath("input.phil").write_text(phil_input)

    result = subprocess.run(
        [shutil.which("dials.combine_experiments"), "input.phil"],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    # load results
    exp = ExperimentListFactory.from_json_file(
        tmp_path / "combined.expt", check_format=False
    )
    ref = flex.reflection_table.from_file(tmp_path / "combined.refl")

    # test the experiments
    assert len(exp) == 103
    assert len(exp.crystals()) == 103
    assert len(exp.beams()) == 1
    assert len(exp.scans()) == 1
    assert len(exp.detectors()) == 1
    assert len(exp.goniometers()) == 1
    for e in exp:
        assert e.imageset is not None

    # test the reflections
    assert len(ref) == 11689

    result = subprocess.run(
        [shutil.which("dials.split_experiments"), "combined.expt", "combined.refl"],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    for i, e in enumerate(exp):
        assert tmp_path.joinpath(f"split_{i:03d}.expt").is_file()
        assert tmp_path.joinpath(f"split_{i:03d}.refl").is_file()

        exp_single = ExperimentListFactory.from_json_file(
            tmp_path / f"split_{i:03d}.expt", check_format=False
        )
        ref_single = flex.reflection_table.from_file(tmp_path / f"split_{i:03d}.refl")

        assert len(exp_single) == 1
        assert exp_single[0].crystal == e.crystal
        assert exp_single[0].beam == e.beam
        assert exp_single[0].detector == e.detector
        assert exp_single[0].scan == e.scan
        assert exp_single[0].goniometer == e.goniometer
        assert exp_single[0].imageset == e.imageset
        assert len(ref_single) == len(ref.select(ref["id"] == i))
        assert ref_single["id"].all_eq(0)

    result = subprocess.run(
        [
            shutil.which("dials.split_experiments"),
            "combined.expt",
            "output.experiments_prefix=test",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    for i in range(len(exp)):
        assert tmp_path.joinpath(f"test_{i:03d}.expt").is_file()

    # Modify a copy of the detector
    detector = copy.deepcopy(exp.detectors()[0])
    panel = detector[0]
    x, y, z = panel.get_origin()
    panel.set_frame(panel.get_fast_axis(), panel.get_slow_axis(), (x, y, z + 10))
    # Set half of the experiments to the new detector
    for i in range(len(exp) // 2):
        exp[i].detector = detector
    exp.as_json(tmp_path / "modded.expt")

    result = subprocess.run(
        [
            shutil.which("dials.split_experiments"),
            "modded.expt",
            "combined.refl",
            "output.experiments_prefix=test_by_detector",
            "output.reflections_prefix=test_by_detector",
            "by_detector=True",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    for i in range(2):
        assert tmp_path.joinpath(f"test_by_detector_{i:03d}.expt").is_file()
        assert tmp_path.joinpath(f"test_by_detector_{i:03d}.refl").is_file()
    assert not tmp_path.joinpath("test_by_detector_002.expt").is_file()
    assert not tmp_path.joinpath("test_by_detector_002.refl").is_file()

    # Now do test when input has identifiers set
    reflections = flex.reflection_table().from_file(tmp_path / "combined.refl")
    explist = ExperimentListFactory.from_json_file(
        tmp_path / "combined.expt", check_format=False
    )
    # set string identifiers as nonconsecutive 0,2,4,6....
    for i, exp in enumerate(explist):
        assert i in reflections["id"]
        reflections.experiment_identifiers()[i] = str(i * 2)
        exp.identifier = str(i * 2)
    reflections.as_file(tmp_path / "assigned.refl")
    explist.as_json(tmp_path / "assigned.expt")

    result = subprocess.run(
        [shutil.which("dials.split_experiments"), "assigned.expt", "assigned.refl"],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    for i in range(len(explist)):
        assert tmp_path.joinpath(f"split_{i:03d}.expt").is_file()
        assert tmp_path.joinpath(f"split_{i:03d}.refl").is_file()

        exp_single = ExperimentListFactory.from_json_file(
            tmp_path / f"split_{i:03d}.expt", check_format=False
        )
        ref_single = flex.reflection_table.from_file(tmp_path / f"split_{i:03d}.refl")

        assert len(exp_single) == 1
        # resets all ids to 0, but keeps mapping to unique identifier.
        # doesn't have to be set to 0 but doing this to keep more consistent with
        # other dials programs
        assert ref_single["id"].all_eq(0)
        assert ref_single.experiment_identifiers()[0] == str(i * 2)

    # update modded experiments to have same identifiers as assigned_experiments
    moddedlist = ExperimentListFactory.from_json_file(
        tmp_path / "modded.expt", check_format=False
    )
    for i, exp in enumerate(moddedlist):
        exp.identifier = str(i * 2)
    moddedlist.as_json(tmp_path / "modded.expt")

    result = subprocess.run(
        [
            shutil.which("dials.split_experiments"),
            "modded.expt",
            "assigned.refl",
            "output.experiments_prefix=test_by_detector",
            "output.reflections_prefix=test_by_detector",
            "by_detector=True",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    # Expect each datasets to have ids from 0..50 with experiment identifiers
    # all kept from before 0,2,4,6,...
    current_exp_id = 0
    for i in range(2):
        assert tmp_path.joinpath(f"test_by_detector_{i:03d}.expt").is_file()
        assert tmp_path.joinpath(f"test_by_detector_{i:03d}.refl").is_file()
        explist = ExperimentListFactory.from_json_file(
            tmp_path / f"test_by_detector_{i:03d}.expt", check_format=False
        )
        refl = flex.reflection_table.from_file(
            tmp_path / f"test_by_detector_{i:03d}.refl"
        )

        for k in range(len(explist)):
            assert refl.experiment_identifiers()[k] == str(current_exp_id)
            current_exp_id += 2

    assert not tmp_path.joinpath("test_by_detector_002.expt").is_file()
    assert not tmp_path.joinpath("test_by_detector_002.refl").is_file()


@pytest.mark.parametrize(
    "with_identifiers,with_reflections",
    [("True", "True"), ("False", "True"), ("True", "False")],
)
def test_combine_clustering(dials_data, tmp_path, with_identifiers, with_reflections):
    """Test with the clustering.use=True option.

    Need to use an integrated dataset for this option.
    """
    data_dir = dials_data("multi_crystal_proteinase_k", pathlib=True)

    input_range = [2, 3, 4, 5, 10]
    if with_identifiers:
        for n, i in enumerate(input_range):
            command = [
                shutil.which("dials.assign_experiment_identifiers"),
                data_dir / f"experiments_{i}.json",
                data_dir / f"reflections_{i}.pickle",
                f"output.experiments={n}.expt",
                f"output.reflections={n}.refl",
            ]
            subprocess.run(command, cwd=tmp_path, capture_output=True)

        phil_input = "\n".join(
            (
                f"  input.experiments={tmp_path / f'{i}.expt'}\n"
                + f"  input.reflections={tmp_path / f'{i}.refl'}"
                if with_reflections
                else ""
            )
            for i in [0, 1, 2, 3, 4]
        )

    else:
        phil_input = "\n".join(
            f"  input.experiments={data_dir}/experiments_{i}.json\n"
            + f"  input.reflections={data_dir}/reflections_{i}.pickle"
            if with_reflections
            else ""
        )

    tmp_path.joinpath("input.phil").write_text(phil_input)

    result = subprocess.run(
        [
            shutil.which("dials.combine_experiments"),
            tmp_path / "input.phil",
            "clustering.use=True",
            "threshold=5",
            "max_clusters=2",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    # this should create two clusters:
    #   combined_cluster_1 (2 expts)
    #   combined_cluster_2 (3 expts)

    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("combined_cluster2.expt").is_file()
    assert tmp_path.joinpath("combined_cluster1.expt").is_file()
    if with_reflections:
        assert tmp_path.joinpath("combined_cluster2.refl").is_file()
        assert tmp_path.joinpath("combined_cluster1.refl").is_file()

    exps = load.experiment_list(tmp_path / "combined_cluster1.expt", check_format=False)
    assert len(exps) == 2
    exps = load.experiment_list(tmp_path / "combined_cluster2.expt", check_format=False)
    assert len(exps) == 3
    if with_reflections:
        refls = flex.reflection_table.from_file(tmp_path / "combined_cluster1.refl")
        assert list(set(refls["id"])) == [0, 1]
        refls = flex.reflection_table.from_file(tmp_path / "combined_cluster2.refl")
        assert list(set(refls["id"])) == [0, 1, 2]


@pytest.fixture
def narrow_wedge_input_with_identifiers(dials_data, tmp_path):
    """Make a fixture to avoid multiple runs of assign identifiers."""
    data_dir = dials_data("polyhedra_narrow_wedges", pathlib=True)
    input_range = [9, 11, 12, 31]
    for n, i in enumerate(input_range):
        command = [
            shutil.which("dials.assign_experiment_identifiers"),
            data_dir / ("sweep_%03d_experiments.json" % i),
            data_dir / ("sweep_%03d_reflections.pickle" % i),
            f"output.experiments={n}.expt",
            f"output.reflections={n}.refl",
        ]
        subprocess.run(command, cwd=tmp_path, capture_output=True)

    phil_input = "\n".join(
        (
            "  input.experiments=%s\n" % (tmp_path / f"{i}.expt")
            + "  input.reflections=%s" % (tmp_path / f"{i}.refl")
        )
        for i, _ in enumerate(input_range)
    )
    return phil_input


@pytest.mark.parametrize("min_refl", ["None", "100"])
@pytest.mark.parametrize("max_refl", ["None", "150"])
def test_min_max_reflections_per_experiment(dials_data, tmp_path, min_refl, max_refl):
    expected_results = {
        ("None", "None"): 10,
        ("None", "150"): 9,
        ("100", "None"): 6,
        ("100", "150"): 5,
    }

    data_dir = dials_data("refinement_test_data", pathlib=True)
    input_phil = (
        f" input.experiments={data_dir}/multi_stills_combined.json\n"
        + f" input.reflections={data_dir}/multi_stills_combined.pickle\n"
        + f" output.min_reflections_per_experiment={min_refl}\n"
        + f" output.max_reflections_per_experiment={max_refl}\n"
    ).format(data_dir, min_refl, max_refl)
    tmp_path.joinpath("input.phil").write_text(input_phil)

    result = subprocess.run(
        [shutil.which("dials.combine_experiments"), "input.phil"],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    # load results
    exp = ExperimentListFactory.from_json_file(
        tmp_path / "combined.expt", check_format=False
    )

    assert len(exp) == expected_results[(min_refl, max_refl)]


@pytest.mark.parametrize("with_identifiers", ["True", "False"])
@pytest.mark.parametrize("method", ["random", "n_refl", "significance_filter"])
def test_combine_nsubset(
    dials_data,
    tmp_path,
    with_identifiers,
    method,
    narrow_wedge_input_with_identifiers,
):
    """Test with the n_subset option."""

    if with_identifiers:
        phil_input = narrow_wedge_input_with_identifiers
    else:
        data_dir = dials_data("polyhedra_narrow_wedges", pathlib=True)
        input_range = [9, 11, 12, 31]
        phil_input = "\n".join(
            (
                "  input.experiments={0}/sweep_%03d_experiments.json\n"
                + "  input.reflections={0}/sweep_%03d_reflections.pickle"
            )
            % (i, i)
            for i in input_range
        ).format(data_dir)

    (tmp_path / "input.phil").write_text(phil_input)

    result = subprocess.run(
        [
            shutil.which("dials.combine_experiments"),
            tmp_path / "input.phil",
            "n_subset=3",
            f"n_subset_method={method}",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "combined.refl").is_file()
    assert (tmp_path / "combined.expt").is_file()

    exps = load.experiment_list(tmp_path / "combined.expt", check_format=False)
    assert len(exps) == 3
    refls = flex.reflection_table.from_file(tmp_path / "combined.refl")
    # Check that order are the same to ensure consistent for historical
    # use of ordered ids to match across datastructures
    assert list(exps.identifiers()) == list(refls.experiment_identifiers().values())
    assert len(set(refls["id"])) == 3
    assert list(set(refls["id"])) == [0, 1, 2]


def test_failed_tolerance_error(dials_data, monkeypatch):
    """Test that we get a sensible error message on tolerance failures"""
    # Select some experiments to use for combining
    data_dir = dials_data("polyhedra_narrow_wedges", pathlib=True)
    jsons = os.path.join(
        data_dir,
        "sweep_{:03d}_{}",
    )
    list_of_elists = [
        load.experiment_list(jsons.format(2, "experiments.json"), check_format=False),
        load.experiment_list(jsons.format(3, "experiments.json"), check_format=False),
    ]
    list_of_elists[0][0].identifier = "0"
    list_of_elists[1][0].identifier = "1"
    list_of_tables = [
        flex.reflection_table.from_file(jsons.format(2, "reflections.pickle")),
        flex.reflection_table.from_file(jsons.format(4, "reflections.pickle")),
    ]

    params = phil_scope.extract()
    params.reference_from_experiment.beam = 0

    # Validate that these pass
    combine_experiments(params, list_of_elists, list_of_tables)

    # Now, alter the beam to check it doesn't pass
    expt2 = list_of_elists[1][0]
    expt2.beam.set_wavelength(expt2.beam.get_wavelength() * 2)

    with pytest.raises(SystemExit) as exc:
        combine_experiments(params, list_of_elists, list_of_tables)
    assert "Beam" in str(exc.value)
    print("Got (expected) error message:", exc.value)


def test_combine_imagesets(dials_data, tmp_path):
    data = dials_data("l_cysteine_dials_output", pathlib=True)
    list_of_elists = [
        load.experiment_list(f, check_format=False)
        for f in sorted(data.glob("*_integrated_experiments.json"))
    ]
    list_of_tables = [
        flex.reflection_table.from_file(f)
        for f in sorted(data.glob("*_integrated.pickle"))
    ]
    params = phil_scope.extract()
    expts, refls = combine_experiments(params, list_of_elists, list_of_tables)
    assert set(refls["imageset_id"]) == {0, 1, 2, 3}

    # Check that we have preserved unindexed reflections for all of these
    assert set(refls.select(refls["id"] == -1)["imageset_id"]) == {0, 1, 2, 3}

    # test combining without reflections
    expts2 = combine_experiments_no_reflections(params, list_of_elists)
    assert len(expts2) == 4
    assert expts2.identifiers() == expts.identifiers()
