"""
Test the command line script dials.scale, for successful completion.
"""

from __future__ import absolute_import, division, print_function

import json
import os
import pytest
import procrunner

from libtbx import phil
from dxtbx.serialize import load
from dxtbx.model.experiment_list import ExperimentList
from dxtbx.model import Crystal, Scan, Beam, Goniometer, Detector, Experiment
from dials.util import Sorry
from dials.array_family import flex
from dials.util.options import OptionParser
from dials.command_line.scale import Script


def run_one_scaling(working_directory, argument_list):
    """Run the dials.scale algorithm."""
    command = ["dials.scale"] + argument_list
    print(command)
    result = procrunner.run(command, working_directory=working_directory)
    assert not result.returncode and not result.stderr
    assert working_directory.join("scaled.expt").check()
    assert working_directory.join("scaled.refl").check()
    assert working_directory.join("scaling.html").check()

    table = flex.reflection_table.from_pickle(
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
    import iotbx.merging_statistics

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


def return_first_arg_side_effect(*args):
    """Side effect for overriding the call to reject_outliers."""
    return args[0]


def test_scale_script_prepare_input():
    """Test prepare_input method of scaling script."""

    # test the components of the scaling script directly with a test reflection
    # table, experiments list and params.

    params, exp, reflections = generate_test_input()
    # try to pass in unequal number of reflections and experiments
    reflections.append(generate_test_reflections())
    with pytest.raises(Sorry):
        _ = Script(params, exp, reflections)

    params, exp, reflections = generate_test_input()
    # Try to use use_datasets when not identifiers set
    params.dataset_selection.use_datasets = ["0"]
    with pytest.raises(Sorry):
        _ = Script(params, exp, reflections)
    # Try to use use_datasets when not identifiers set
    params.dataset_selection.use_datasets = None
    params.dataset_selection.exclude_datasets = ["0"]
    with pytest.raises(Sorry):
        _ = Script(params, exp, reflections)

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
    params.dataset_selection.use_datasets = ["0"]
    params, exp, script_reflections = Script.prepare_input(params, exp, reflections)

    assert len(script_reflections) == 1

    # Try again, this time excluding
    params, exp, reflections = generate_test_input(n=2)
    exp[0].identifier = "0"
    reflections[0].experiment_identifiers()[0] = "0"
    exp[1].identifier = "1"
    reflections[1].experiment_identifiers()[0] = "1"
    params.dataset_selection.exclude_datasets = ["0"]
    params, exp, script_reflections = Script.prepare_input(params, exp, reflections)

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
    with pytest.raises(Sorry):
        _ = Script.prepare_input(params, exp, reflections)

    # Test cutting data
    params, exp, reflections = generate_test_input(n=1)
    params.cut_data.d_min = 1.5
    params, _, script_reflections = Script.prepare_input(params, exp, reflections)
    r = script_reflections[0]
    assert list(r.get_flags(r.flags.user_excluded_in_scaling)) == [
        False,
        False,
        True,
        True,
    ]
    params.cut_data.d_max = 2.25
    params, _, script_reflections = Script.prepare_input(params, exp, reflections)
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
    _, __, script_reflections = Script.prepare_input(params, exp, reflections)
    r = script_reflections[0]
    assert list(r.get_flags(r.flags.user_excluded_in_scaling)) == [
        True,
        False,
        False,
        False,
    ]


def test_scale_physical(dials_regression, tmpdir):
    """Test standard scaling of one dataset."""

    data_dir = os.path.join(dials_regression, "xia2-28")
    refl_1 = os.path.join(data_dir, "20_integrated.pickle")
    expt_1 = os.path.join(data_dir, "20_integrated_experiments.json")
    extra_args = [
        "model=physical",
        "merged_mtz=merged.mtz",
        "optimise_errors=False",
        "intensity_choice=profile",
        "unmerged_mtz=unmerged.mtz",
        "use_free_set=1",
        "outlier_rejection=simple",
    ]
    run_one_scaling(tmpdir, [refl_1, expt_1] + extra_args)
    assert tmpdir.join("unmerged.mtz").check()
    assert tmpdir.join("merged.mtz").check()

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
        result.overall.r_pim < 0.0255
    )  # at 14/08/18, value was 0.023, at 07/02/19 was 0.0243
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
        "concurrent=False",
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
    # test the 'stats_only' option
    extra_args = ["stats_only=True"]
    run_one_scaling(tmpdir, ["scaled.refl", "scaled.expt"] + extra_args)
    # test the 'export_mtz_only' option
    extra_args = [
        "export_mtz_only=True",
        "unmerged_mtz=test_1.mtz",
        "merged_mtz=test_2.mtz",
    ]
    run_one_scaling(tmpdir, ["scaled.refl", "scaled.expt"] + extra_args)
    assert tmpdir.join("test_1.mtz").check()
    assert tmpdir.join("test_2.mtz").check()


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
        "optimise_errors=False",
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
    assert result.overall.r_pim < 0.17  # 17/05/19 was 0.1525
    assert result.overall.cc_one_half > 0.95  # 17/05/19 was 0.9722, 29/07/19 was 0.9557
    assert result.overall.n_obs > 51400  # 17/05/19 was 51560, 29/07/19 was 51493
    # for this dataset, expect to have two regions excluded - last 5 images of
    # datasets _4 & _5
    with open(tmpdir.join("analysis_results.json").strpath) as f:
        analysis_results = json.load(f)
    assert analysis_results["cycle_results"]["1"]["image_ranges_removed"] == [
        [[21, 25], 4]
    ]
    assert analysis_results["cycle_results"]["2"]["image_ranges_removed"] == [
        [[21, 25], 3],
        [[21, 25], 5],
    ]
    assert analysis_results["termination_reason"] == "no_more_removed"


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
        "optimise_errors=False",
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
    assert analysis_results["cycle_results"]["1"]["removed_datasets"] == ["4"]


def test_scale_optimise_errors(dials_regression, tmpdir):
    """Test standard scaling of one dataset with error optimisation."""
    data_dir = os.path.join(dials_regression, "xia2-28")
    refl = os.path.join(data_dir, "20_integrated.pickle")
    expt = os.path.join(data_dir, "20_integrated_experiments.json")
    extra_args = ["model=physical", "optimise_errors=True"]
    run_one_scaling(tmpdir, [refl, expt] + extra_args)


def test_scale_array(dials_regression, tmpdir):
    """Test a standard dataset - ideally needs a large dataset or full matrix
    round may fail. Currently turning off absorption term to avoid
    overparameterisation and failure of full matrix minimisation."""

    data_dir = os.path.join(dials_regression, "xia2-28")
    refl = os.path.join(data_dir, "20_integrated.pickle")
    expt = os.path.join(data_dir, "20_integrated_experiments.json")
    extra_args = ["model=array", "absorption_term=0", "full_matrix=0"]

    run_one_scaling(tmpdir, [refl, expt] + extra_args)


def test_multi_scale(dials_regression, tmpdir):
    """Test standard scaling of two datasets."""

    data_dir = os.path.join(dials_regression, "xia2-28")
    refl_1 = os.path.join(data_dir, "20_integrated.pickle")
    expt_1 = os.path.join(data_dir, "20_integrated_experiments.json")
    refl_2 = os.path.join(data_dir, "25_integrated.pickle")
    expt_2 = os.path.join(data_dir, "25_integrated_experiments.json")
    extra_args = [
        "unmerged_mtz=unmerged.mtz",
        "optimise_errors=False",
        "intensity_choice=profile",
        "outlier_rejection=simple",
    ]

    run_one_scaling(tmpdir, [refl_1, refl_2, expt_1, expt_2] + extra_args)

    # Now inspect output, check it hasn't changed drastically, or if so verify
    # that the new behaviour is more correct and update test accordingly.
    result = get_merging_stats(tmpdir.join("unmerged.mtz").strpath)
    expected_nobs = 5460
    assert abs(result.overall.n_obs - expected_nobs) < 30
    assert result.overall.r_pim < 0.0221  # at 22/10/18, value was 0.22037
    assert result.overall.cc_one_half > 0.9975  # at 07/08/18, value was 0.99810
    print(result.overall.r_pim)
    print(result.overall.cc_one_half)

    # run again, optimising errors, and continuing from where last run left off.
    extra_args = [
        "optimise_errors=True",
        "unmerged_mtz=unmerged.mtz",
        "check_consistent_indexing=True",
    ]
    run_one_scaling(tmpdir, ["scaled.refl", "scaled.expt"] + extra_args)
    # Now inspect output, check it hasn't changed drastically, or if so verify
    # that the new behaviour is more correct and update test accordingly.
    # Note: error optimisation currently appears to give worse results here!
    result = get_merging_stats(tmpdir.join("unmerged.mtz").strpath)
    expected_nobs = 5411
    print(result.overall.r_pim)
    print(result.overall.cc_one_half)
    assert abs(result.overall.n_obs - expected_nobs) < 100
    assert result.overall.r_pim < 0.023  # at 07/08/18, value was 0.022722
    assert result.overall.cc_one_half > 0.9965  # at 07/08/18, value was 0.996925


def test_multi_scale_exclude_images(dials_regression, tmpdir):
    """Test scaling of multiple dataset with image exclusion."""
    data_dir = os.path.join(dials_regression, "xia2-28")
    refl_1 = os.path.join(data_dir, "20_integrated.pickle")
    expt_1 = os.path.join(data_dir, "20_integrated_experiments.json")
    refl_2 = os.path.join(data_dir, "25_integrated.pickle")
    expt_2 = os.path.join(data_dir, "25_integrated_experiments.json")
    # Expect this dataset to be given batches 1-1800 and 1901-3600
    # Try excluding last two hundred batches
    extra_args = [
        "optimise_errors=False",
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
        "optimise_errors=False",
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


def test_targeted_scaling(dials_regression, tmpdir, dials_data):
    """Test the targeted scaling workflow."""
    location = dials_data("l_cysteine_4_sweeps_scaled")
    target_refl = location.join("scaled_35.refl").strpath
    target_expt = location.join("scaled_35.expt").strpath

    data_dir = os.path.join(dials_regression, "xia2-28")
    refl_1 = os.path.join(data_dir, "20_integrated.pickle")
    expt_1 = os.path.join(data_dir, "20_integrated_experiments.json")
    refl_2 = os.path.join(data_dir, "25_integrated.pickle")
    expt_2 = os.path.join(data_dir, "25_integrated_experiments.json")

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


def test_incremental_scale_workflow(dials_regression, tmpdir):
    """Try scale, cosym, scale, cosym, scale."""
    data_dir = os.path.join(dials_regression, "xia2-28")
    refl_1 = os.path.join(data_dir, "20_integrated.pickle")
    expt_1 = os.path.join(data_dir, "20_integrated_experiments.json")

    run_one_scaling(tmpdir, [refl_1, expt_1])

    # test order also - first new file before scaled
    refl_2 = os.path.join(data_dir, "25_integrated.pickle")
    expt_2 = os.path.join(data_dir, "25_integrated_experiments.json")
    command = ["dials.cosym", refl_2, expt_2, "scaled.refl", "scaled.expt"]

    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("symmetrized.expt").check()
    assert tmpdir.join("symmetrized.refl").check()

    # now try scaling again to check everything is okay
    run_one_scaling(tmpdir, ["symmetrized.refl", "symmetrized.expt"])

    # test order also - first scaled file then new file
    refl_2 = os.path.join(data_dir, "30_integrated.pickle")
    expt_2 = os.path.join(data_dir, "30_integrated_experiments.json")
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
        ("multi", "absorption_term", None),
        ("multi", "model", "physical array"),
    ],
)
def test_scale_cross_validate(
    dials_regression, tmpdir, mode, parameter, parameter_values
):
    """Test standard scaling of one dataset."""
    data_dir = os.path.join(dials_regression, "xia2-28")
    refl = os.path.join(data_dir, "20_integrated.pickle")
    expt = os.path.join(data_dir, "20_integrated_experiments.json")
    extra_args = [
        "cross_validation_mode=%s" % mode,
        "nfolds=2",
        "full_matrix=0",
        "optimise_errors=0",
    ]
    if parameter:
        extra_args += ["parameter=%s" % parameter]
    if parameter_values:
        extra_args += ["parameter_values=%s" % parameter_values]
    command = ["dials.scale", refl, expt] + extra_args
    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
