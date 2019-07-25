from __future__ import absolute_import, division, print_function

import mock
import pytest
from dials.util.options import OptionParser
from dials.algorithms.scaling.cross_validation.crossvalidator import (
    DialsScaleCrossValidator,
)
from dials.algorithms.scaling.cross_validation.cross_validate import cross_validate
from libtbx import phil


def generated_param():
    """Generate a param phil scope."""
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


def test_crossvalidator():
    """Use the dials.scale cross validator to test the general methods of the
    parent class."""
    experiments = []
    reflections = []

    # test the create_results_dict method
    crossvalidator = DialsScaleCrossValidator(experiments, reflections)
    crossvalidator.create_results_dict(1)
    assert len(crossvalidator.results_dict) == 1
    assert (
        len(crossvalidator.results_dict[0])
        == len(crossvalidator.results_metadata["names"]) + 1
    )
    crossvalidator = DialsScaleCrossValidator(experiments, reflections)
    crossvalidator.create_results_dict(2)
    assert len(crossvalidator.results_dict) == 2
    assert (
        len(crossvalidator.results_dict[0])
        == len(crossvalidator.results_metadata["names"]) + 1
    )
    assert (
        len(crossvalidator.results_dict[1])
        == len(crossvalidator.results_metadata["names"]) + 1
    )

    # now test the add configuration
    keys = ("model",)
    values = (["a", "b"],)
    crossvalidator.set_results_dict_configuration(keys, values)
    assert crossvalidator.results_dict[0]["configuration"][0] == "model=a"
    assert crossvalidator.results_dict[1]["configuration"][0] == "model=b"

    # now try adding results to dict
    results_1 = [1.0, 1.5, 0.5, 3.0, 4.0, 1.0, 0.1, 0.2, 0.1]
    results_2 = [1.0, 2.0, 1.0, 3.0, 2.0, 2.0, 0.3, 0.4, 0.1]
    results_3 = [1.0, 2.0, 1.0, 3.0, 3.0, 0.0, 0.1, 0.3, 0.2]
    results_4 = [1.0, 2.0, 1.0, 3.0, 3.5, 0.5, 0.3, 0.5, 0.2]
    crossvalidator.add_results_to_results_dict(0, results_1)
    assert crossvalidator.results_dict[0]["work Rpim"] == [1.0]
    assert crossvalidator.results_dict[0]["free Rpim"] == [1.5]
    assert crossvalidator.results_dict[0]["work CC1/2"] == [3.0]
    assert crossvalidator.results_dict[0]["free CC1/2"] == [4.0]
    # add some more so that can test interpreting results
    crossvalidator.add_results_to_results_dict(0, results_2)
    crossvalidator.add_results_to_results_dict(1, results_3)
    crossvalidator.add_results_to_results_dict(1, results_4)

    # Now try interpreting results - check that values are calculated correctly
    st = crossvalidator.interpret_results()
    r1 = [
        "model=a",
        "mean",
        "1.0",
        "1.75*",
        "0.75*",
        "3.0",
        "3.0",
        "1.5",
        "0.2",
        "0.3",
        "0.1*",
    ]
    r2 = [
        "",
        "std dev",
        "0.0",
        "0.25",
        "0.25",
        "0.0",
        "1.0",
        "0.5",
        "0.1",
        "0.1",
        "0.0",
    ]
    r3 = [
        "model=b",
        "mean",
        "1.0",
        "2.0",
        "1.0",
        "3.0",
        "3.25*",
        "0.25*",
        "0.2",
        "0.4*",
        "0.2",
    ]
    r4 = [
        "",
        "std dev",
        "0.0",
        "0.0",
        "0.0",
        "0.0",
        "0.25",
        "0.25",
        "0.1",
        "0.1",
        "0.0",
    ]
    assert st._rows[0] == r1
    assert st._rows[1] == r2
    assert st._rows[2] == r3
    assert st._rows[3] == r4


def test_dialsscalecrossvalidator():
    """Test the methods of the dials.scale cross validator"""
    experiments = []
    reflections = []

    def mock_script():
        script = mock.MagicMock()
        script.scaler.final_rmsds = [1.0, 2.0, 3.0, 4.0]
        return script

    # test get results from script
    crossvalidator = DialsScaleCrossValidator(experiments, reflections)
    script = mock_script()
    results = crossvalidator.get_results_from_script(script)
    assert results == [1.0, 2.0, 3.0, 4.0]

    params = generated_param()
    params.scaling_options.free_set_percentage = 20.0
    # test get free set offset
    fsp = crossvalidator.get_free_set_percentage(params)
    assert fsp == 20.0

    # test set free set offset
    params = crossvalidator.set_free_set_offset(params, 5)
    assert params.scaling_options.free_set_offset == 5

    # test get/set parameters
    assert crossvalidator.get_parameter_type("model") == "choice"
    assert crossvalidator.get_parameter_type("absorption_term") == "bool"
    assert crossvalidator.get_parameter_type("lmax") == "int"
    assert crossvalidator.get_parameter_type("decay_interval") == "float"

    params = crossvalidator.set_parameter(params, "model", "KB")
    assert params.model == "KB"
    params = crossvalidator.set_parameter(params, "decay_interval", 50.0)
    assert params.parameterisation.decay_interval == 50.0
    params = crossvalidator.set_parameter(params, "lmax", 10)
    assert params.parameterisation.lmax == 10
    params = crossvalidator.set_parameter(params, "optimise_errors", False)
    assert params.weighting.optimise_errors is False
    params = crossvalidator.set_parameter(params, "d_min", 1.8)
    assert params.cut_data.d_min == 1.8
    params = crossvalidator.set_parameter(params, "outlier_zmax", 7.53)
    assert params.scaling_options.outlier_zmax == 7.53
    with pytest.raises(AssertionError):
        _ = crossvalidator.set_parameter(params, "bad_parameter", 7.53)

    # defer testing of run_script to command line tests


def test_cross_validate_script():
    """Test the script, mocking the run_script and interpret results calls"""

    param = generated_param()
    crossvalidator = DialsScaleCrossValidator([], [])

    # test expected error raise due to unspecified parameter
    param.cross_validation.cross_validation_mode = "multi"
    with pytest.raises(ValueError):
        cross_validate(param, crossvalidator)

    # test single mode
    param.cross_validation.cross_validation_mode = "single"
    param.cross_validation.nfolds = 2
    fpath = "dials.algorithms.scaling.cross_validation."
    with mock.patch(
        fpath + "crossvalidator.DialsScaleCrossValidator.run_script"
    ) as mock_run_script:
        with mock.patch(
            fpath + "crossvalidator.DialsScaleCrossValidator.interpret_results"
        ) as mock_interpret:
            cross_validate(param, crossvalidator)
            assert mock_run_script.call_count == 2
            assert mock_interpret.call_count == 1

    # test multi mode
    param = generated_param()
    param.cross_validation.cross_validation_mode = "multi"
    param.cross_validation.nfolds = 2
    fpath = "dials.algorithms.scaling.cross_validation."
    with mock.patch(
        fpath + "crossvalidator.DialsScaleCrossValidator.run_script"
    ) as mock_run_script:
        with mock.patch(
            fpath + "crossvalidator.DialsScaleCrossValidator.interpret_results"
        ) as mock_interpret:
            param.cross_validation.parameter = "absorption_term"
            cross_validate(param, crossvalidator)
            assert mock_run_script.call_count == 4
            assert mock_interpret.call_count == 1

            param.cross_validation.parameter = "decay_interval"
            with pytest.raises(ValueError):
                cross_validate(param, crossvalidator)

            param.cross_validation.parameter = "absorption_term"
            param.cross_validation.parameter_values = ["True", "False"]
            cross_validate(param, crossvalidator)
            assert mock_run_script.call_count == 8
            assert mock_interpret.call_count == 2

            param.cross_validation.parameter = "decay_interval"
            param.cross_validation.parameter_values = ["5.0", "10.0"]
            cross_validate(param, crossvalidator)
            assert mock_run_script.call_count == 12
            assert mock_interpret.call_count == 3

            param.cross_validation.parameter = "model"
            param.cross_validation.parameter_values = ["array", "physical"]
            cross_validate(param, crossvalidator)
            assert mock_run_script.call_count == 16
            assert mock_interpret.call_count == 4

            param.cross_validation.parameter = "lmax"
            param.cross_validation.parameter_values = ["4", "6"]
            cross_validate(param, crossvalidator)
            assert mock_run_script.call_count == 20
            assert mock_interpret.call_count == 5

            param.cross_validation.parameter = "bad_interval"
            with pytest.raises(ValueError):
                cross_validate(param, crossvalidator)

            param.cross_validation.cross_validation_mode = "bad"
            with pytest.raises(ValueError):
                cross_validate(param, crossvalidator)
