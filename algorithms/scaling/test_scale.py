# coding: utf-8

"""
Test the command line script dials.scale, for successful completion.
"""

from __future__ import absolute_import, division, print_function

import json

import iotbx.merging_statistics
import pytest
import procrunner
from cctbx import uctbx
import iotbx.mtz
from libtbx import phil
from dxtbx.serialize import load
from dxtbx.model.experiment_list import ExperimentList
from dxtbx.model import Crystal, Scan, Beam, Goniometer, Detector, Experiment
from dials.array_family import flex
from dials.util.options import OptionParser
from dials.algorithms.scaling.algorithm import ScalingAlgorithm, prepare_input
from dials.command_line import merge, report, scale


def run_one_scaling(working_directory, argument_list):
    """Run the dials.scale algorithm."""
    command = ["dials.scale"] + argument_list
    print(command)
    result = procrunner.run(command, working_directory=working_directory)
    print(result.stderr)
    assert not result.returncode and not result.stderr
    assert working_directory.join("scaled.expt").check()
    assert working_directory.join("scaled.refl").check()
    assert working_directory.join("dials.scale.html").check()

    table = flex.reflection_table.from_file(
        working_directory.join("scaled.refl").strpath
    )

    assert "inverse_scale_factor" in table
    assert "inverse_scale_factor_variance" in table


def get_merging_stats(
    scaled_unmerged_mtz,
    anomalous=False,
    n_bins=20,
    use_internal_variance=False,
    eliminate_sys_absent=False,
    data_labels=None,
):
    """Return a merging statistics result from an mtz file."""

    i_obs = iotbx.merging_statistics.select_data(
        scaled_unmerged_mtz, data_labels=data_labels
    )
    i_obs = i_obs.customized_copy(anomalous_flag=False, info=i_obs.info())
    result = iotbx.merging_statistics.dataset_statistics(
        i_obs=i_obs,
        n_bins=n_bins,
        anomalous=anomalous,
        use_internal_variance=use_internal_variance,
        eliminate_sys_absent=eliminate_sys_absent,
    )
    return result


def generated_exp(n=1):
    """Generate an experiment list with two experiments."""
    experiments = ExperimentList()
    exp_dict = {
        "__id__": "crystal",
        "real_space_a": [1.0, 0.0, 0.0],
        "real_space_b": [0.0, 1.0, 0.0],
        "real_space_c": [0.0, 0.0, 2.0],
        "space_group_hall_symbol": " C 2y",
    }
    crystal = Crystal.from_dict(exp_dict)
    scan = Scan(image_range=[0, 90], oscillation=[0.0, 1.0])
    beam = Beam(s0=(0.0, 0.0, 1.01))
    goniometer = Goniometer((1.0, 0.0, 0.0))
    detector = Detector()
    experiments.append(
        Experiment(
            beam=beam,
            scan=scan,
            goniometer=goniometer,
            detector=detector,
            crystal=crystal,
        )
    )
    if n > 1:
        for _ in range(n - 1):
            experiments.append(
                Experiment(
                    beam=beam,
                    scan=scan,
                    goniometer=goniometer,
                    detector=detector,
                    crystal=crystal,
                )
            )
    return experiments


def generated_param():
    """Generate the default scaling parameters object."""
    phil_scope = phil.parse(
        """
      include scope dials.command_line.scale.phil_scope
  """,
        process_includes=True,
    )

    optionparser = OptionParser(phil=phil_scope, check_format=False)
    parameters, _ = optionparser.parse_args(
        args=[], quick_parse=True, show_diff_phil=False
    )
    return parameters


def generate_test_reflections():
    """Generate a small reflection table."""
    reflections = flex.reflection_table()
    reflections["intensity.sum.value"] = flex.double([1.0, 2.0, 3.0, 4.0])
    reflections["intensity.sum.variance"] = flex.double([1.0, 2.0, 3.0, 4.0])
    reflections["miller_index"] = flex.miller_index(
        [(0, 0, 1), (0, 0, 1), (0, 0, 2), (0, 0, 2)]
    )
    reflections["xyzobs.px.value"] = flex.vec3_double(
        [(0, 0, 1), (0, 0, 2), (0, 0, 3), (0, 0, 4)]
    )
    reflections["id"] = flex.int([0, 0, 0, 0])
    return reflections


def generate_test_input(n=1):
    """Generate params, exps and refls."""
    reflections = []
    for _ in range(n):
        reflections.append(generate_test_reflections())
    return generated_param(), generated_exp(n), reflections


def test_scale_script_prepare_input():
    """Test prepare_input method of scaling script."""

    # test the components of the scaling script directly with a test reflection
    # table, experiments list and params.

    params, exp, reflections = generate_test_input()
    # try to pass in unequal number of reflections and experiments
    reflections.append(generate_test_reflections())
    with pytest.raises(ValueError):
        _ = ScalingAlgorithm(params, exp, reflections)

    params, exp, reflections = generate_test_input()
    # Try to use use_datasets when not identifiers set
    params.dataset_selection.use_datasets = [0]
    with pytest.raises(ValueError):
        _ = ScalingAlgorithm(params, exp, reflections)
    # Try to use use_datasets when not identifiers set
    params.dataset_selection.use_datasets = None
    params.dataset_selection.exclude_datasets = [0]
    with pytest.raises(ValueError):
        _ = ScalingAlgorithm(params, exp, reflections)

    # Now make two experiments with identifiers and select on them
    params, exp, reflections = generate_test_input(n=2)
    exp[0].identifier = "0"
    reflections[0].experiment_identifiers()[0] = "0"
    exp[1].identifier = "1"
    reflections[1].experiment_identifiers()[0] = "1"
    list1 = ExperimentList().append(exp[0])
    list2 = ExperimentList().append(exp[1])
    reflections[0].assert_experiment_identifiers_are_consistent(list1)
    reflections[1].assert_experiment_identifiers_are_consistent(list2)
    params.dataset_selection.use_datasets = [0]
    params, exp, script_reflections = prepare_input(params, exp, reflections)

    assert len(script_reflections) == 1

    # Try again, this time excluding
    params, exp, reflections = generate_test_input(n=2)
    exp[0].identifier = "0"
    reflections[0].experiment_identifiers()[0] = "0"
    exp[1].identifier = "1"
    reflections[1].experiment_identifiers()[1] = "1"
    params.dataset_selection.exclude_datasets = [0]
    params, exp, script_reflections = prepare_input(params, exp, reflections)

    assert len(script_reflections) == 1
    assert script_reflections[0] is reflections[1]

    # Try having two unequal space groups
    params, exp, reflections = generate_test_input(n=2)
    exp_dict = {
        "__id__": "crystal",
        "real_space_a": [1.0, 0.0, 0.0],
        "real_space_b": [0.0, 1.0, 0.0],
        "real_space_c": [0.0, 0.0, 2.0],
        "space_group_hall_symbol": " P 1",
    }
    crystal = Crystal.from_dict(exp_dict)
    exp[0].crystal = crystal
    with pytest.raises(ValueError):
        _ = prepare_input(params, exp, reflections)

    # Test cutting data
    params, exp, reflections = generate_test_input(n=1)
    params.cut_data.d_min = 1.5
    params, _, script_reflections = prepare_input(params, exp, reflections)
    r = script_reflections[0]
    assert list(r.get_flags(r.flags.user_excluded_in_scaling)) == [
        False,
        False,
        True,
        True,
    ]
    params.cut_data.d_max = 2.25
    params, _, script_reflections = prepare_input(params, exp, reflections)
    r = script_reflections[0]
    assert list(r.get_flags(r.flags.user_excluded_in_scaling)) == [
        False,
        False,
        True,
        True,
    ]

    params, exp, reflections = generate_test_input(n=1)
    reflections[0]["partiality"] = flex.double([0.5, 0.8, 1.0, 1.0])
    params.cut_data.partiality_cutoff = 0.75
    _, __, script_reflections = prepare_input(params, exp, reflections)
    r = script_reflections[0]
    assert list(r.get_flags(r.flags.user_excluded_in_scaling)) == [
        True,
        False,
        False,
        False,
    ]


def test_targeted_scaling_against_mtz(dials_data, tmpdir):
    """Test targeted scaling against an mtz generated with dials.scale."""
    location = dials_data("l_cysteine_4_sweeps_scaled")
    refl = location.join("scaled_35.refl").strpath
    expt = location.join("scaled_35.expt").strpath
    command = ["dials.scale", refl, expt, "unmerged_mtz=unmerged.mtz"]

    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("scaled.expt").check()
    assert tmpdir.join("scaled.refl").check()
    assert tmpdir.join("unmerged.mtz").check()

    refl = location.join("scaled_30.refl").strpath
    expt = location.join("scaled_30.expt").strpath
    target_mtz = tmpdir.join("unmerged.mtz").strpath
    command = ["dials.scale", refl, expt, "target_mtz=%s" % target_mtz]

    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("scaled.expt").check()
    assert tmpdir.join("scaled.refl").check()
    expts = load.experiment_list(tmpdir.join("scaled.expt").strpath, check_format=False)
    assert len(expts) == 1


@pytest.mark.parametrize(
    "option",
    [
        None,
        "reflection_selection.method=random",
        "reflection_selection.method=intensity_ranges",
        "reflection_selection.method=use_all",
        "intensity_choice=sum",
        "intensity_choice=profile",
    ],
)
def test_scale_single_dataset_with_options(dials_data, tmpdir, option):
    """Test different non-default command-line options with a single dataset."""
    data_dir = dials_data("l_cysteine_dials_output")
    refl_1 = data_dir / "20_integrated.pickle"
    expt_1 = data_dir / "20_integrated_experiments.json"
    args = [refl_1, expt_1]
    if option:
        args.append(option)
    run_one_scaling(tmpdir, args)


@pytest.fixture
def vmxi_protk_reindexed(dials_data, tmpdir):
    """Reindex the protk data to be in the correct space group."""
    location = dials_data("vmxi_proteinase_k_sweeps")

    command = [
        "dials.reindex",
        location.join("experiments_0.json"),
        location.join("reflections_0.pickle"),
        "space_group=P422",
    ]
    procrunner.run(command, working_directory=tmpdir)
    return tmpdir.join("reindexed.expt"), tmpdir.join("reindexed.refl")


@pytest.mark.parametrize(
    ("options", "expected", "tolerances"),
    [
        (["error_model=None"], None, None),
        (
            ["error_model=basic", "basic.minimisation=individual"],
            (0.73711, 0.04720),
            (0.05, 0.005),
        ),
        (["error_model.basic.a=0.73711"], (0.73711, 0.04720), (1e-6, 0.005)),
        (["error_model.basic.b=0.04720"], (0.73711, 0.04720), (0.05, 1e-6)),
        (
            ["error_model.basic.b=0.02", "error_model.basic.a=1.5"],
            (1.50, 0.02),
            (1e-6, 1e-6),
        ),
        (
            ["error_model=basic", "basic.minimisation=regression"],
            (0.995, 0.051),
            (0.05, 0.005),
        ),
        (
            ["error_model.basic.a=0.99", "basic.minimisation=regression"],
            (0.99, 0.051),
            (1e-6, 0.005),
        ),
        (
            ["error_model.basic.b=0.051", "basic.minimisation=regression"],
            (0.99, 0.051),
            (0.05, 1e-6),
        ),
    ],
)
def test_error_model_options(
    vmxi_protk_reindexed, tmpdir, options, expected, tolerances
):
    """Test different non-default command-line options with a single dataset.

    Current values taken at 14.11.19"""
    expt_1, refl_1 = vmxi_protk_reindexed
    args = [refl_1, expt_1] + [o for o in options]
    run_one_scaling(tmpdir, args)
    expts = load.experiment_list(tmpdir.join("scaled.expt").strpath, check_format=False)
    config = expts[0].scaling_model.configdict
    if not expected:
        assert "error_model_parameters" not in config
    else:
        params = expts[0].scaling_model.configdict["error_model_parameters"]
        print(list(params))
        assert params[0] == pytest.approx(expected[0], abs=tolerances[0])
        assert params[1] == pytest.approx(expected[1], abs=tolerances[1])


@pytest.mark.parametrize(
    "option",
    [
        None,
        "combine.joint_analysis=False",
        "reflection_selection.method=quasi_random",
        "reflection_selection.method=random",
        "reflection_selection.method=intensity_ranges",
        "reflection_selection.method=use_all",
    ],
)
def test_scale_multiple_datasets_with_options(dials_data, tmpdir, option):
    """Test different non-defaul command-line options with multiple datasets."""
    data_dir = dials_data("l_cysteine_dials_output")
    refl_1 = data_dir / "20_integrated.pickle"
    expt_1 = data_dir / "20_integrated_experiments.json"
    refl_2 = data_dir / "25_integrated.pickle"
    expt_2 = data_dir / "25_integrated_experiments.json"
    args = [refl_1, expt_1, refl_2, expt_2]
    if option:
        args.append(option)
    run_one_scaling(tmpdir, args)


def test_scale_physical(dials_data, tmpdir):
    """Test standard scaling of one dataset."""

    data_dir = dials_data("l_cysteine_dials_output")
    refl_1 = data_dir / "20_integrated.pickle"
    expt_1 = data_dir / "20_integrated_experiments.json"
    extra_args = [
        "model=physical",
        "merged_mtz=merged.mtz",
        "error_model=None",
        "intensity_choice=profile",
        "unmerged_mtz=unmerged.mtz",
        "crystal_name=foo",
        "project_name=bar",
        "use_free_set=1",
        "outlier_rejection=simple",
        "json=scaling.json",
    ]
    run_one_scaling(tmpdir, [refl_1, expt_1] + extra_args)
    assert tmpdir.join("unmerged.mtz").check()
    assert tmpdir.join("merged.mtz").check()
    assert tmpdir.join("scaling.json").check()
    for f in ("unmerged.mtz", "merged.mtz"):
        mtz_obj = iotbx.mtz.object(tmpdir.join(f).strpath)
        assert mtz_obj.crystals()[1].name() == "foo"
        assert mtz_obj.crystals()[1].project_name() == "bar"

    # Now inspect output, check it hasn't changed drastically, or if so verify
    # that the new behaviour is more correct and update test accordingly.
    result = get_merging_stats(tmpdir.join("unmerged.mtz").strpath)
    print(result.overall.r_pim, result.overall.cc_one_half, result.overall.n_obs)
    assert result.overall.r_pim < 0.0255  # at 30/01/19, value was 0.02410
    assert result.overall.cc_one_half > 0.9955  # at 30/01/19, value was 0.9960
    assert result.overall.n_obs > 2300  # at 30/01/19, was 2320

    # Try running again with the merged.mtz as a target, to trigger the
    # target_mtz option
    extra_args.append("target_mtz=merged.mtz")
    run_one_scaling(tmpdir, [refl_1, expt_1] + extra_args)
    result = get_merging_stats(tmpdir.join("unmerged.mtz").strpath)
    assert (
        result.overall.r_pim < 0.026
    )  # at 17/07/19 was 0.256 after updates to merged mtz export
    assert (
        result.overall.cc_one_half > 0.9955
    )  # at 14/08/18, value was 0.999, at 07/02/19 was 0.9961
    assert result.overall.n_obs > 2300  # at 07/01/19, was 2321, at 07/02/19 was 2321

    # run again with the concurrent scaling option turned off and the 'standard'
    # outlier rejection
    extra_args = [
        "model=physical",
        "merged_mtz=merged.mtz",
        "unmerged_mtz=unmerged.mtz",
        "use_free_set=1",
        "outlier_rejection=standard",
        "refinement_order=consecutive",
        "intensity_choice=combine",
    ]
    run_one_scaling(tmpdir, [refl_1, expt_1] + extra_args)

    # Now inspect output, check it hasn't changed drastically, or if so verify
    # that the new behaviour is more correct and update test accordingly.
    result = get_merging_stats(tmpdir.join("unmerged.mtz").strpath)
    assert (
        result.overall.r_pim < 0.024
    )  # at 07/01/19, value was 0.02372, at 30/01/19 was 0.021498
    assert (
        result.overall.cc_one_half > 0.995
    )  # at 07/01/19, value was 0.99568, at 30/01/19 was 0.9961
    assert result.overall.n_obs > 2300  # at 07/01/19, was 2336, at 22/05/19 was 2311


def test_scale_and_filter_image_group_mode(dials_data, tmpdir):
    """Test the scale and filter command line program."""
    location = dials_data("multi_crystal_proteinase_k")

    command = [
        "dials.scale",
        "filtering.method=deltacchalf",
        "stdcutoff=3.0",
        "mode=image_group",
        "max_cycles=6",
        "d_min=1.4",
        "group_size=5",
        "unmerged_mtz=unmerged.mtz",
        "scale_and_filter_results=analysis_results.json",
        "error_model=None",
    ]
    for i in [1, 2, 3, 4, 5, 7, 10]:
        command.append(location.join("experiments_" + str(i) + ".json").strpath)
        command.append(location.join("reflections_" + str(i) + ".pickle").strpath)

    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("scaled.refl").check()
    assert tmpdir.join("scaled.expt").check()
    assert tmpdir.join("analysis_results.json").check()
    result = get_merging_stats(tmpdir.join("unmerged.mtz").strpath)
    assert result.overall.r_pim < 0.17  # 03/02/20 was 0.160
    assert result.overall.cc_one_half > 0.95  # 03/02/20 was 0.961
    assert result.overall.n_obs > 50000  # 03/02/20 was 50213

    with open(tmpdir.join("analysis_results.json").strpath) as f:
        analysis_results = json.load(f)
    assert analysis_results["cycle_results"]["1"]["image_ranges_removed"] == [
        [[16, 24], 4]
    ]
    assert analysis_results["cycle_results"]["2"]["image_ranges_removed"] == [
        [[17, 24], 3]
    ]
    assert analysis_results["cycle_results"]["3"]["image_ranges_removed"] == [
        [[21, 25], 5]
    ]
    assert analysis_results["termination_reason"] == "max_percent_removed"


def test_scale_and_filter_image_group_single_dataset(dials_data, tmpdir):
    """Test the scale and filter deltacchalf.mode=image_group on a
       single data set."""
    data_dir = dials_data("l_cysteine_dials_output")
    command = [
        "dials.scale",
        data_dir / "20_integrated.pickle",
        data_dir / "20_integrated_experiments.json",
        "filtering.method=deltacchalf",
        "stdcutoff=3.0",
        "mode=image_group",
        "max_cycles=1",
        "scale_and_filter_results=analysis_results.json",
        "error_model=None",
    ]
    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("scaled.refl").check()
    assert tmpdir.join("scaled.expt").check()
    assert tmpdir.join("analysis_results.json").check()

    with open(tmpdir.join("analysis_results.json").strpath) as f:
        analysis_results = json.load(f)
    assert analysis_results["cycle_results"]["1"]["image_ranges_removed"] == []
    assert len(analysis_results["cycle_results"].keys()) == 1
    assert analysis_results["termination_reason"] == "no_more_removed"


def test_scale_dose_decay_model(dials_data, tmpdir):
    """Test the scale and filter command line program."""
    location = dials_data("multi_crystal_proteinase_k")
    command = ["dials.scale", "d_min=2.0", "model=dose_decay"]
    for i in [1, 2, 3, 4, 5, 7, 10]:
        command.append(location.join("experiments_" + str(i) + ".json").strpath)
        command.append(location.join("reflections_" + str(i) + ".pickle").strpath)

    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("scaled.refl").check()
    assert tmpdir.join("scaled.expt").check()
    assert tmpdir.join("dials.scale.html").check()
    expts = load.experiment_list(tmpdir.join("scaled.expt").strpath, check_format=False)
    assert expts[0].scaling_model.id_ == "dose_decay"


def test_scale_best_unit_cell_d_min(dials_data, tmpdir):
    location = dials_data("multi_crystal_proteinase_k")
    best_unit_cell = uctbx.unit_cell((42, 42, 39, 90, 90, 90))
    d_min = 2
    command = [
        "dials.scale",
        "best_unit_cell=%g,%g,%g,%g,%g,%g" % best_unit_cell.parameters(),
        "d_min=%g" % d_min,
        "unmerged_mtz=unmerged.mtz",
    ]
    for i in [1, 2, 3, 4, 5, 7, 10]:
        command.append(location.join("experiments_" + str(i) + ".json").strpath)
        command.append(location.join("reflections_" + str(i) + ".pickle").strpath)
    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("scaled.refl").check()
    assert tmpdir.join("scaled.expt").check()
    stats = get_merging_stats(tmpdir.join("unmerged.mtz").strpath)
    assert stats.overall.d_min >= d_min
    assert stats.crystal_symmetry.unit_cell().parameters() == pytest.approx(
        best_unit_cell.parameters()
    )


def test_scale_and_filter_dataset_mode(dials_data, tmpdir):
    """Test the scale and filter command line program."""
    location = dials_data("multi_crystal_proteinase_k")
    command = [
        "dials.scale",
        "filtering.method=deltacchalf",
        "stdcutoff=1.0",
        "mode=dataset",
        "max_cycles=2",
        "d_min=1.4",
        "output.reflections=filtered.refl",
        "scale_and_filter_results=analysis_results.json",
        "unmerged_mtz=unmerged.mtz",
        "error_model=None",
    ]
    for i in [1, 2, 3, 4, 5, 7, 10]:
        command.append(location.join("experiments_" + str(i) + ".json").strpath)
        command.append(location.join("reflections_" + str(i) + ".pickle").strpath)

    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("filtered.refl").check()
    assert tmpdir.join("scaled.expt").check()
    assert tmpdir.join("analysis_results.json").check()
    with open(tmpdir.join("analysis_results.json").strpath) as f:
        analysis_results = json.load(f)

    assert analysis_results["cycle_results"]["1"]["removed_datasets"] == [
        analysis_results["initial_expids_and_image_ranges"][4][0]
    ]


def test_scale_array(dials_data, tmpdir):
    """Test a standard dataset - ideally needs a large dataset or full matrix
    round may fail. Currently turning off absorption term to avoid
    overparameterisation and failure of full matrix minimisation."""

    data_dir = dials_data("l_cysteine_dials_output")
    refl = data_dir / "20_integrated.pickle"
    expt = data_dir / "20_integrated_experiments.json"
    extra_args = ["model=array", "array.absorption_correction=0", "full_matrix=0"]

    run_one_scaling(tmpdir, [refl, expt] + extra_args)


def test_multi_scale(dials_data, tmpdir):
    """Test standard scaling of two datasets."""

    data_dir = dials_data("l_cysteine_dials_output")
    refl_1 = data_dir / "20_integrated.pickle"
    expt_1 = data_dir / "20_integrated_experiments.json"
    refl_2 = data_dir / "25_integrated.pickle"
    expt_2 = data_dir / "25_integrated_experiments.json"
    extra_args = [
        "unmerged_mtz=unmerged.mtz",
        "error_model=None",
        "intensity_choice=profile",
        "outlier_rejection=simple",
    ]

    run_one_scaling(tmpdir, [refl_1, refl_2, expt_1, expt_2] + extra_args)

    # Now inspect output, check it hasn't changed drastically, or if so verify
    # that the new behaviour is more correct and update test accordingly.
    result = get_merging_stats(tmpdir.join("unmerged.mtz").strpath)
    expected_nobs = 5526  # 19/06/20
    assert abs(result.overall.n_obs - expected_nobs) < 30
    assert result.overall.r_pim < 0.0221  # at 22/10/18, value was 0.22037
    assert result.overall.cc_one_half > 0.9975  # at 07/08/18, value was 0.99810
    print(result.overall.r_pim)
    print(result.overall.cc_one_half)

    # run again, optimising errors, and continuing from where last run left off.
    extra_args = [
        "error_model=basic",
        "unmerged_mtz=unmerged.mtz",
        "check_consistent_indexing=True",
    ]
    run_one_scaling(tmpdir, ["scaled.refl", "scaled.expt"] + extra_args)
    # Now inspect output, check it hasn't changed drastically, or if so verify
    # that the new behaviour is more correct and update test accordingly.
    # Note: error optimisation currently appears to give worse results here!
    result = get_merging_stats(tmpdir.join("unmerged.mtz").strpath)
    expected_nobs = 5560  # 22/06/20
    print(result.overall.r_pim)
    print(result.overall.cc_one_half)
    assert abs(result.overall.n_obs - expected_nobs) < 100
    assert result.overall.r_pim < 0.016  # at #22/06/20, value was 0.015
    assert result.overall.cc_one_half > 0.997  # at #22/06/20, value was 0.999


def test_multi_scale_exclude_images(dials_data, tmpdir):
    """Test scaling of multiple dataset with image exclusion."""
    data_dir = dials_data("l_cysteine_dials_output")
    refl_1 = data_dir / "20_integrated.pickle"
    expt_1 = data_dir / "20_integrated_experiments.json"
    refl_2 = data_dir / "25_integrated.pickle"
    expt_2 = data_dir / "25_integrated_experiments.json"
    # Expect this dataset to be given batches 1-1800 and 1901-3600
    # Try excluding last two hundred batches
    extra_args = [
        "error_model=None",
        "intensity_choice=profile",
        "outlier_rejection=simple",
        "exclude_images=0:1601:1800",
        "exclude_images=1:1501:1700",
    ]

    run_one_scaling(tmpdir, [refl_1, refl_2, expt_1, expt_2] + extra_args)

    scaling_models = load.experiment_list(
        tmpdir.join("scaled.expt").strpath, check_format=False
    ).scaling_models()
    assert scaling_models[0].configdict["valid_image_range"] == [1, 1600]
    assert scaling_models[1].configdict["valid_image_range"] == [1, 1500]
    assert pytest.approx(scaling_models[0].configdict["valid_osc_range"], [0, 160.0])
    assert pytest.approx(scaling_models[1].configdict["valid_osc_range"], [-145.0, 5.0])

    # Run again, excluding some more from one run.
    extra_args = [
        "error_model=None",
        "intensity_choice=profile",
        "outlier_rejection=simple",
        "exclude_images=0:1401:1600",
    ]
    run_one_scaling(tmpdir, ["scaled.refl", "scaled.expt"] + extra_args)
    scaling_models = load.experiment_list(
        tmpdir.join("scaled.expt").strpath, check_format=False
    ).scaling_models()
    assert scaling_models[0].configdict["valid_image_range"] == [1, 1400]
    assert scaling_models[1].configdict["valid_image_range"] == [1, 1500]
    assert pytest.approx(scaling_models[0].configdict["valid_osc_range"], [0, 140.0])
    assert pytest.approx(scaling_models[1].configdict["valid_osc_range"], [-145.0, 5.0])


def test_targeted_scaling(dials_data, tmpdir):
    """Test the targeted scaling workflow."""
    location = dials_data("l_cysteine_4_sweeps_scaled")
    target_refl = location.join("scaled_35.refl").strpath
    target_expt = location.join("scaled_35.expt").strpath

    data_dir = dials_data("l_cysteine_dials_output")
    refl_1 = data_dir / "20_integrated.pickle"
    expt_1 = data_dir / "20_integrated_experiments.json"
    refl_2 = data_dir / "25_integrated.pickle"
    expt_2 = data_dir / "25_integrated_experiments.json"

    # Do targeted scaling, use this as a chance to test the KB model as well.
    extra_args = ["model=KB"]
    run_one_scaling(tmpdir, [target_refl, refl_1, target_expt, expt_1] + extra_args)

    scaled_exp = tmpdir.join("scaled.expt").strpath
    scaled_refl = tmpdir.join("scaled.refl").strpath
    experiments_list = load.experiment_list(scaled_exp, check_format=False)
    assert len(experiments_list.scaling_models()) == 2
    assert experiments_list.scaling_models()[0].id_ == "physical"
    assert experiments_list.scaling_models()[1].id_ == "KB"

    extra_args = ["model=KB", "only_target=True"]
    run_one_scaling(tmpdir, [refl_2, scaled_refl, expt_2, scaled_exp] + extra_args)


def test_incremental_scale_workflow(dials_data, tmpdir):
    """Try scale, cosym, scale, cosym, scale."""
    data_dir = dials_data("l_cysteine_dials_output")
    refl_1 = data_dir / "20_integrated.pickle"
    expt_1 = data_dir / "20_integrated_experiments.json"

    run_one_scaling(tmpdir, [refl_1, expt_1])

    # test order also - first new file before scaled
    refl_2 = data_dir / "25_integrated.pickle"
    expt_2 = data_dir / "25_integrated_experiments.json"
    command = ["dials.cosym", refl_2, expt_2, "scaled.refl", "scaled.expt"]

    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("symmetrized.expt").check()
    assert tmpdir.join("symmetrized.refl").check()

    # now try scaling again to check everything is okay
    run_one_scaling(tmpdir, ["symmetrized.refl", "symmetrized.expt"])

    # test order also - first scaled file then new file
    refl_2 = data_dir / "30_integrated.pickle"
    expt_2 = data_dir / "30_integrated_experiments.json"
    command = [
        "dials.cosym",
        "scaled.refl",
        "scaled.expt",
        refl_1,
        expt_2,
        "output.reflections=symmetrized.refl",
        "output.experiments=symmetrized.expt",
    ]
    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("symmetrized.expt").check()
    assert tmpdir.join("symmetrized.refl").check()

    # now try scaling again to check everything is okay
    args = ["dials.scale", "symmetrized.refl", "symmetrized.expt"]
    command = " ".join(args)
    run_one_scaling(tmpdir, ["symmetrized.refl", "symmetrized.expt"])


@pytest.mark.parametrize(
    ("mode", "parameter", "parameter_values"),
    [
        ("single", None, None),
        ("multi", "physical.absorption_correction", None),
        ("multi", "model", "physical array"),
    ],
)
def test_scale_cross_validate(dials_data, tmpdir, mode, parameter, parameter_values):
    """Test standard scaling of one dataset."""
    data_dir = dials_data("l_cysteine_dials_output")
    refl = data_dir / "20_integrated.pickle"
    expt = data_dir / "20_integrated_experiments.json"
    extra_args = [
        "cross_validation_mode=%s" % mode,
        "nfolds=2",
        "full_matrix=0",
        "error_model=None",
    ]
    if parameter:
        extra_args += ["parameter=%s" % parameter]
    if parameter_values:
        extra_args += ["parameter_values=%s" % parameter_values]
    command = ["dials.scale", refl, expt] + extra_args
    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr


@pytest.mark.xfail(
    reason="test state leakage, cf. https://github.com/dials/dials/issues/1271",
)
def test_few_reflections(dials_data):
    u"""
    Test that dials.symmetry does something sensible if given few reflections.

    Use some example integrated data generated from two ten-image 1° sweeps.  These
    each contain a few dozen integrated reflections.

    Also test the behaviour of dials.merge and dials.report on the output.

    By suppressing the output from dials.scale and dials.report, we obviate the need to
    run in a temporary directory.

    Args:
        dials_data: DIALS custom Pytest fixture for access to test data.
    """
    # Get the input experiment lists and reflection tables.
    data_dir = dials_data("l_cysteine_dials_output")
    experiments = ExperimentList.from_file(data_dir / "11_integrated.expt")
    experiments.extend(ExperimentList.from_file(data_dir / "23_integrated.expt"))
    refls = "11_integrated.refl", "23_integrated.refl"
    reflections = [flex.reflection_table.from_file(data_dir / refl) for refl in refls]

    # Get and use the default parameters for dials.scale, suppressing HTML output.
    scale_params = scale.phil_scope.fetch(
        source=phil.parse("output.html=None")
    ).extract()
    # Does what it says on the tin.  Run scaling.
    scaled_expt, scaled_refl = scale.run_scaling(scale_params, experiments, reflections)

    # Get and use the default parameters for dials.merge.
    merge_params = merge.phil_scope.fetch(source=phil.parse("")).extract()
    # Run dials.merge on the scaling output.
    merge.merge_data_to_mtz(merge_params, scaled_expt, [scaled_refl])

    # Get and use the default parameters for dials.report, suppressing HTML output.
    report_params = report.phil_scope.fetch(
        source=phil.parse("output.html=None")
    ).extract()
    # Get an Analyser object, which does the dials.report stuff.
    analyse = report.Analyser(
        report_params,
        grid_size=report_params.grid_size,
        centroid_diff_max=report_params.centroid_diff_max,
    )
    # Run dials.report on scaling output.
    analyse(scaled_refl, scaled_expt)
