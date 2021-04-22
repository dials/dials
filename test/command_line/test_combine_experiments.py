"""
Test combination of multiple experiments and reflections files.
"""

from __future__ import absolute_import, division, print_function

import copy
import os

import procrunner
import pytest

from dxtbx.model.experiment_list import ExperimentListFactory
from dxtbx.serialize import load

import dials.command_line.combine_experiments as combine_experiments
from dials.array_family import flex


def test(dials_regression, run_in_tmpdir):
    data_dir = os.path.join(
        dials_regression, "refinement_test_data", "multi_narrow_wedges"
    )

    input_range = list(range(2, 49))
    for i in (8, 10, 15, 16, 34, 39, 45):
        input_range.remove(i)

    phil_input = "\n".join(
        (
            "  input.experiments={0}/data/sweep_%03d/experiments.json\n"
            + "  input.reflections={0}/data/sweep_%03d/reflections.pickle"
        )
        % (i, i)
        for i in input_range
    )

    # assert phil_input == "\n" + phil_input2 + "\n "

    input_phil = (
        phil_input.format(data_dir)
        + """
 reference_from_experiment.beam=0
 reference_from_experiment.scan=0
 reference_from_experiment.goniometer=0
 reference_from_experiment.detector=0
 """
    )

    with open("input.phil", "w") as phil_file:
        phil_file.writelines(input_phil)

    result = procrunner.run(["dials.combine_experiments", "input.phil"])
    assert not result.returncode and not result.stderr

    # load results
    exp = ExperimentListFactory.from_json_file("combined.expt", check_format=False)
    ref = flex.reflection_table.from_file("combined.refl")

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

    result = procrunner.run(
        ["dials.split_experiments", "combined.expt", "combined.refl"]
    )
    assert not result.returncode and not result.stderr

    for i, e in enumerate(exp):
        assert os.path.exists("split_%03d.expt" % i)
        assert os.path.exists("split_%03d.refl" % i)

        exp_single = ExperimentListFactory.from_json_file(
            "split_%03d.expt" % i, check_format=False
        )
        ref_single = flex.reflection_table.from_file("split_%03d.refl" % i)

        assert len(exp_single) == 1
        assert exp_single[0].crystal == e.crystal
        assert exp_single[0].beam == e.beam
        assert exp_single[0].detector == e.detector
        assert exp_single[0].scan == e.scan
        assert exp_single[0].goniometer == e.goniometer
        assert exp_single[0].imageset == e.imageset
        assert len(ref_single) == len(ref.select(ref["id"] == i))
        assert ref_single["id"].all_eq(0)

    result = procrunner.run(
        ["dials.split_experiments", "combined.expt", "output.experiments_prefix=test"]
    )
    assert not result.returncode and not result.stderr

    for i in range(len(exp)):
        assert os.path.exists("test_%03d.expt" % i)

    # Modify a copy of the detector
    detector = copy.deepcopy(exp.detectors()[0])
    panel = detector[0]
    x, y, z = panel.get_origin()
    panel.set_frame(panel.get_fast_axis(), panel.get_slow_axis(), (x, y, z + 10))
    # Set half of the experiments to the new detector
    for i in range(len(exp) // 2):
        exp[i].detector = detector
    exp.as_json("modded.expt")

    result = procrunner.run(
        [
            "dials.split_experiments",
            "modded.expt",
            "combined.refl",
            "output.experiments_prefix=test_by_detector",
            "output.reflections_prefix=test_by_detector",
            "by_detector=True",
        ]
    )
    assert not result.returncode and not result.stderr

    for i in range(2):
        assert os.path.exists("test_by_detector_%03d.expt" % i)
        assert os.path.exists("test_by_detector_%03d.refl" % i)
    assert not os.path.exists("test_by_detector_%03d.expt" % 2)
    assert not os.path.exists("test_by_detector_%03d.refl" % 2)

    # Now do test when input has identifiers set
    reflections = flex.reflection_table().from_file("combined.refl")
    explist = ExperimentListFactory.from_json_file("combined.expt", check_format=False)
    # set string identifiers as nonconsecutive 0,2,4,6....
    for i, exp in enumerate(explist):
        assert i in reflections["id"]
        reflections.experiment_identifiers()[i] = str(i * 2)
        exp.identifier = str(i * 2)
    reflections.as_file("assigned.refl")
    explist.as_json("assigned.expt")

    result = procrunner.run(
        ["dials.split_experiments", "assigned.expt", "assigned.refl"]
    )
    assert not result.returncode and not result.stderr

    for i in range(len(explist)):
        assert os.path.exists("split_%03d.expt" % i)
        assert os.path.exists("split_%03d.refl" % i)

        exp_single = ExperimentListFactory.from_json_file(
            "split_%03d.expt" % i, check_format=False
        )
        ref_single = flex.reflection_table.from_file("split_%03d.refl" % i)

        assert len(exp_single) == 1
        # resets all ids to 0, but keeps mapping to unique identifier.
        # doesn't have to be set to 0 but doing this to keep more consistent with
        # other dials programs
        assert ref_single["id"].all_eq(0)
        assert ref_single.experiment_identifiers()[0] == str(i * 2)

    # update modded experiments to have same identifiers as assigned_experiments
    moddedlist = ExperimentListFactory.from_json_file("modded.expt", check_format=False)
    for i, exp in enumerate(moddedlist):
        exp.identifier = str(i * 2)
    moddedlist.as_json("modded.expt")

    result = procrunner.run(
        [
            "dials.split_experiments",
            "modded.expt",
            "assigned.refl",
            "output.experiments_prefix=test_by_detector",
            "output.reflections_prefix=test_by_detector",
            "by_detector=True",
        ]
    )
    assert not result.returncode and not result.stderr

    # Expect each datasets to have ids from 0..50 with experiment identifiers
    # all kept from before 0,2,4,6,...
    current_exp_id = 0
    for i in range(2):
        assert os.path.exists("test_by_detector_%03d.expt" % i)
        assert os.path.exists("test_by_detector_%03d.refl" % i)
        explist = ExperimentListFactory.from_json_file(
            "test_by_detector_%03d.expt" % i, check_format=False
        )
        refl = flex.reflection_table.from_file("test_by_detector_%03d.refl" % i)

        for k in range(len(explist)):
            assert refl.experiment_identifiers()[k] == str(current_exp_id)
            current_exp_id += 2

    assert not os.path.exists("test_by_detector_%03d.expt" % 2)
    assert not os.path.exists("test_by_detector_%03d.refl" % 2)


@pytest.mark.parametrize("with_identifiers", ["True", "False"])
def test_combine_clustering(dials_data, tmpdir, with_identifiers):
    """Test with the clustering.use=True option.

    Need to use an integrated dataset for this option.
    """
    data_dir = dials_data("multi_crystal_proteinase_k")

    input_range = [2, 3, 4, 5, 10]
    if with_identifiers:
        for n, i in enumerate(input_range):
            command = [
                "dials.assign_experiment_identifiers",
                data_dir.strpath + "/experiments_%s.json" % i,
                data_dir.strpath + "/reflections_%s.pickle" % i,
                "output.experiments=%s.expt" % n,
                "output.reflections=%s.refl" % n,
            ]
            procrunner.run(command, working_directory=tmpdir)

        phil_input = "\n".join(
            (
                "  input.experiments=%s\n" % tmpdir.join("%s.expt" % i)
                + "  input.reflections=%s" % tmpdir.join("%s.refl" % i)
            )
            for i in [0, 1, 2, 3, 4]
        )
    else:
        phil_input = "\n".join(
            (
                "  input.experiments={0}/experiments_%s.json\n"
                + "  input.reflections={0}/reflections_%s.pickle"
            )
            % (i, i)
            for i in input_range
        ).format(data_dir.strpath)

    with open(tmpdir.join("input.phil").strpath, "w") as phil_file:
        phil_file.writelines(phil_input)

    result = procrunner.run(
        [
            "dials.combine_experiments",
            tmpdir.join("input.phil").strpath,
            "clustering.use=True",
            "threshold=5",
            "max_clusters=2",
        ],
        working_directory=tmpdir,
    )
    # this should create two clusters:
    #   combined_cluster_1 (2 expts)
    #   combined_cluster_2 (3 expts)

    assert not result.returncode and not result.stderr
    assert tmpdir.join("combined_cluster2.refl").check()
    assert tmpdir.join("combined_cluster2.expt").check()
    assert tmpdir.join("combined_cluster1.refl").check()
    assert tmpdir.join("combined_cluster1.expt").check()

    exps = load.experiment_list(
        tmpdir.join("combined_cluster1.expt").strpath, check_format=False
    )
    assert len(exps) == 2
    refls = flex.reflection_table.from_file(tmpdir.join("combined_cluster1.refl"))
    assert list(set(refls["id"])) == [0, 1]
    exps = load.experiment_list(
        tmpdir.join("combined_cluster2.expt").strpath, check_format=False
    )
    assert len(exps) == 3
    refls = flex.reflection_table.from_file(tmpdir.join("combined_cluster2.refl"))
    assert list(set(refls["id"])) == [0, 1, 2]


@pytest.fixture
def narrow_wedge_input_with_identifiers(dials_regression, tmpdir):
    """Make a fixture to avoid multiple runs of assign identifiers."""
    data_dir = os.path.join(
        dials_regression, "refinement_test_data", "multi_narrow_wedges"
    )
    input_range = [9, 11, 12, 31]
    for n, i in enumerate(input_range):
        command = [
            "dials.assign_experiment_identifiers",
            os.path.join(data_dir, "data/sweep_%03d/experiments.json" % i),
            os.path.join(data_dir, "data/sweep_%03d/reflections.pickle" % i),
            "output.experiments=%s.expt" % n,
            "output.reflections=%s.refl" % n,
        ]
        procrunner.run(command, working_directory=tmpdir)

    phil_input = "\n".join(
        (
            "  input.experiments=%s\n" % tmpdir.join("%s.expt" % i)
            + "  input.reflections=%s" % tmpdir.join("%s.refl" % i)
        )
        for i, _ in enumerate(input_range)
    )
    return phil_input


@pytest.mark.parametrize("min_refl", ["None", "100"])
@pytest.mark.parametrize("max_refl", ["None", "150"])
def test_min_max_reflections_per_experiment(
    dials_regression, run_in_tmpdir, min_refl, max_refl
):

    expected_results = {
        ("None", "None"): 10,
        ("None", "150"): 9,
        ("100", "None"): 6,
        ("100", "150"): 5,
    }

    data_dir = os.path.join(dials_regression, "refinement_test_data", "multi_stills")
    input_phil = (
        " input.experiments={0}/combined_experiments.json\n"
        + " input.reflections={0}/combined_reflections.pickle\n"
        + " output.min_reflections_per_experiment={1}\n"
        + " output.max_reflections_per_experiment={2}\n"
    ).format(data_dir, min_refl, max_refl)
    with open("input.phil", "w") as phil_file:
        phil_file.writelines(input_phil)

    result = procrunner.run(["dials.combine_experiments", "input.phil"])
    assert not result.returncode and not result.stderr

    # load results
    exp = ExperimentListFactory.from_json_file("combined.expt", check_format=False)

    assert len(exp) == expected_results[(min_refl, max_refl)]


@pytest.mark.parametrize("with_identifiers", ["True", "False"])
@pytest.mark.parametrize("method", ["random", "n_refl", "significance_filter"])
def test_combine_nsubset(
    dials_regression,
    tmpdir,
    with_identifiers,
    method,
    narrow_wedge_input_with_identifiers,
):
    """Test with the n_subset option."""

    if with_identifiers:
        phil_input = narrow_wedge_input_with_identifiers
    else:
        data_dir = os.path.join(
            dials_regression, "refinement_test_data", "multi_narrow_wedges"
        )
        input_range = [9, 11, 12, 31]
        phil_input = "\n".join(
            (
                "  input.experiments={0}/data/sweep_%03d/experiments.json\n"
                + "  input.reflections={0}/data/sweep_%03d/reflections.pickle"
            )
            % (i, i)
            for i in input_range
        ).format(data_dir)

    with open(tmpdir.join("input.phil").strpath, "w") as phil_file:
        phil_file.writelines(phil_input)

    result = procrunner.run(
        [
            "dials.combine_experiments",
            tmpdir.join("input.phil").strpath,
            "n_subset=3",
            "n_subset_method=%s" % method,
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("combined.refl").check()
    assert tmpdir.join("combined.expt").check()

    exps = load.experiment_list(
        tmpdir.join("combined.expt").strpath, check_format=False
    )
    assert len(exps) == 3
    refls = flex.reflection_table.from_file(tmpdir.join("combined.refl"))
    # Check that order are the same to ensure consistent for historical
    # use of ordered ids to match across datastructures
    assert list(exps.identifiers()) == list(refls.experiment_identifiers().values())
    assert len(set(refls["id"])) == 3
    assert list(set(refls["id"])) == [0, 1, 2]


def test_failed_tolerance_error(dials_regression, monkeypatch):
    """Test that we get a sensible error message on tolerance failures"""
    # Select some experiments to use for combining
    jsons = os.path.join(
        dials_regression,
        "refinement_test_data",
        "multi_narrow_wedges",
        "data",
        "sweep_{:03d}",
        "{}",
    )
    files = [
        jsons.format(2, "experiments.json"),
        jsons.format(2, "reflections.pickle"),
        jsons.format(3, "experiments.json"),
        jsons.format(3, "reflections.pickle"),
    ]

    # Use the combine script
    script = combine_experiments.Script()
    # Disable writing output
    monkeypatch.setattr(script, "_save_output", lambda *args: None)
    # Parse arguments and configure
    params, options = script.parser.parse_args(files)
    params.reference_from_experiment.beam = 0

    # Validate that these pass
    script.run_with_preparsed(params, options)

    # Now, alter the beam to check it doesn't pass
    exp_2 = params.input.experiments[1].data[0]
    exp_2.beam.set_wavelength(exp_2.beam.get_wavelength() * 2)

    with pytest.raises(SystemExit) as exc:
        script.run_with_preparsed(params, options)
    assert "Beam" in str(exc.value)
    print("Got (expected) error message:", exc.value)


def test_mpi_combine(dials_data, tmpdir):
    pytest.importorskip("mpi4py")
    data_dir = dials_data("multi_crystal_proteinase_k")
    cmd_ref = [
        "dials.combine_experiments",
        data_dir.join("experiments_*.*").strpath,
        data_dir.join("reflection_*.*").strpath,
        "output.experiments=" + tmpdir.join("ref.expt").strpath,
        "output.reflections=" + tmpdir.join("ref.refl").strpath,
    ]
    cmd_test = [
        "mpirun",
        "-n",
        "3",
        "dials.combine_experiments",
        data_dir.join("experiments_*.*").strpath,
        data_dir.join("reflection_*.*").strpath,
        "input.experiments_suffix=.json",
        "input.reflections_suffix=.pickle",
        "output.experiments=" + tmpdir.join("test.expt").strpath,
        "output.reflections=" + tmpdir.join("test.refl").strpath,
    ]

    result_ref = procrunner.run(cmd_ref, working_directory=tmpdir)
    result_test = procrunner.run(cmd_test, working_directory=tmpdir)
    assert not result_ref.returncode and not result_ref.stderr
    assert not result_test.returncode and not result_test.stderr

    exps_ref = load.experiment_list(tmpdir.join("ref.expt").strpath)
    exps_test = load.experiment_list(tmpdir.join("test.expt").strpath)
    assert len(exps_ref) == len(exps_test)

    refls_ref = flex.reflection_table.from_file(tmpdir.join("ref.refl"))
    refls_test = flex.reflection_table.from_file(tmpdir.join("test.refl"))
    assert len(refls_ref) == len(refls_test)
