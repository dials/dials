import os
import mock
from collections import OrderedDict
from dials.array_family import flex
from dials.algorithms.scaling.model.model import KBScalingModel
from dxtbx.model import Experiment
from dials.util.observer import Subject
from dials.algorithms.scaling.observers import (
    register_default_scaling_observers,
    ScalingModelObserver,
    ScalingOutlierObserver,
    ScalingHTMLGenerator,
    ErrorModelObserver,
    MergingStatisticsObserver,
    ScalingSummaryGenerator,
)


def test_register_scaling_observers():
    """Test the registering of the standard scaling observers."""

    class Scaler(Subject):

        """Test scaler class"""

        def __init__(self):
            super(Scaler, self).__init__(
                events=[
                    "performed_scaling",
                    "performed_outlier_rejection",
                    "performed_error_analysis",
                ]
            )

    class TestScript(Subject):

        """Test script class"""

        def __init__(self):
            super(TestScript, self).__init__(
                events=["merging_statistics", "run_script"]
            )
            self.scaler = Scaler()

    script = TestScript()
    register_default_scaling_observers(script)
    assert script.get_observers("merging_statistics") == {
        MergingStatisticsObserver(): MergingStatisticsObserver().update
    }
    assert script.get_observers("run_script") == {
        ScalingHTMLGenerator(): ScalingHTMLGenerator().make_scaling_html,
        ScalingSummaryGenerator(): ScalingSummaryGenerator().print_scaling_summary,
    }
    assert script.scaler.get_observers("performed_scaling") == {
        ScalingModelObserver(): ScalingModelObserver().update
    }
    assert script.scaler.get_observers("performed_outlier_rejection") == {
        ScalingOutlierObserver(): ScalingOutlierObserver().update
    }


def test_ScalingHTMLGenerator(run_in_tmpdir):
    """Test the scaling html generator."""
    script = mock.Mock()
    script.params.output.html = "test.html"

    # Test that ScalingHTMLGenerator works if all data is empty.
    observer = ScalingHTMLGenerator()
    observer.make_scaling_html(script)
    assert os.path.exists("test.html")


def test_ScalingModelObserver():
    """Test that the observer correctly logs data when passed a scaler."""

    KB_dict = {
        "__id__": "KB",
        "is_scaled": True,
        "scale": {
            "n_parameters": 1,
            "parameters": [0.5],
            "est_standard_devs": [0.05],
            "null_parameter_value": 1,
        },
        "configuration_parameters": {"corrections": ["scale"]},
    }

    KBmodel = KBScalingModel.from_dict(KB_dict)
    experiment = Experiment()
    experiment.scaling_model = KBmodel
    experiment.identifier = "0"

    scaler1 = mock.Mock()
    scaler1.experiment = experiment
    scaler1.active_scalers = None

    observer = ScalingModelObserver()
    observer.update(scaler1)
    assert observer.data["0"] == KB_dict

    msg = observer.return_model_error_summary()
    assert msg != ""

    mock_func = mock.Mock()
    mock_func.return_value = {"plot": {}}

    with mock.patch(
        "dials.algorithms.scaling.observers.plot_scaling_models", new=mock_func
    ):
        observer.make_plots()
        assert mock_func.call_count == 1
        assert mock_func.call_args_list == [mock.call(KB_dict)]

    experiment2 = Experiment()
    experiment2.scaling_model = KBmodel
    experiment2.identifier = "1"
    scaler2 = mock.Mock()
    scaler2.experiment = experiment2
    scaler2.active_scalers = None

    multiscaler = mock.Mock()
    multiscaler.active_scalers = [scaler1, scaler2]
    observer.data = {}
    observer.update(multiscaler)
    assert observer.data["0"] == KB_dict
    assert observer.data["1"] == KB_dict

    mock_func = mock.Mock()
    mock_func.return_value = {"plot": {}}

    with mock.patch(
        "dials.algorithms.scaling.observers.plot_scaling_models", new=mock_func
    ):
        r = observer.make_plots()
        assert mock_func.call_count == 2
        assert mock_func.call_args_list == [mock.call(KB_dict), mock.call(KB_dict)]
        assert r == {"scaling_model": {"plot_0": {}, "plot_1": {}}}


def test_ScalingOutlierObserver():
    """Test that the observer correctly logs data when passed a scaler."""

    mock_detector = mock.Mock()
    mock_detector.get_image_size.return_value = [100.0, 200.0]
    mock_scan = mock.Mock()
    mock_scan.get_oscillation.return_value = [0.0, 1.0]
    mock_scan.get_oscillation_range.return_value = [0.0, 4.0]

    scaler = mock.Mock()
    scaler.suitable_refl_for_scaling_sel = flex.bool(4, True)
    scaler.outliers = flex.bool([False, True, False, False])
    scaler.reflection_table = {
        "xyzobs.px.value": flex.vec3_double(
            [(0, 1, 2), (10, 11, 12), (20, 21, 22), (30, 31, 32)]
        )
    }
    scaler.experiment.identifier = "0"
    scaler.experiment.detector = [mock_detector]
    scaler.experiment.scan = mock_scan
    scaler.active_scalers = None

    observer = ScalingOutlierObserver()
    observer.update(scaler)
    assert observer.data == {
        "0": {
            "x": [10],
            "y": [11],
            "z": [12],
            "image_size": [100.0, 200.0],
            "z_range": [0.0, 4.0],
        }
    }

    mock_func = mock.Mock()
    mock_func.return_value = {"outlier_xy_positions": {}, "outliers_vs_z": {}}

    with mock.patch("dials.algorithms.scaling.observers.plot_outliers", new=mock_func):
        r = observer.make_plots()
        assert mock_func.call_count == 1
        assert mock_func.call_args_list == [mock.call(observer.data["0"])]
        assert r == {
            "outlier_plots": OrderedDict(
                [("outlier_plot_0", {}), ("outlier_plot_z0", {})]
            )
        }


def test_ErrorModelObserver():
    """Test that the observer correctly logs data when passed a scaler."""
    delta_hl = flex.double(range(10))

    scaler = mock.Mock()
    scaler.experiment.scaling_model.error_model.delta_hl = delta_hl
    scaler.active_scalers = None

    observer = ErrorModelObserver()
    observer.update(scaler)
    assert observer.data["delta_hl"] == list(delta_hl)
    mock_func = mock.Mock()
    mock_func.return_value = {"plot": {}}

    with mock.patch(
        "dials.algorithms.scaling.observers.normal_probability_plot", new=mock_func
    ):
        r = observer.make_plots()
        assert mock_func.call_count == 1
        assert mock_func.call_args_list == [mock.call(observer.data)]
        assert r == {"normal_prob_plot": {"plot": {}}}


def test_MergingStatisticsObserver():
    """Test that the observer correctly logs data when passed a script."""
    script = mock.Mock()
    script.merging_statistics_result = "result"
    script.scaled_miller_array.space_group.return_value.is_centric.return_value = True
    observer = MergingStatisticsObserver()
    observer.update(script)

    assert observer.data == {"statistics": "result", "is_centric": True}

    mock_func = mock.Mock()
    mock_func.return_value = "return_tables"
    mock_func_2 = mock.Mock()
    mock_func_2.return_value = "cchalfplot"

    with mock.patch(
        "dials.algorithms.scaling.observers.statistics_tables", new=mock_func
    ):
        with mock.patch(
            "dials.algorithms.scaling.observers.cc_one_half_plot", new=mock_func_2
        ):
            r = observer.make_plots()
            assert mock_func.call_count == 1
            assert mock_func.call_args_list == [mock.call(observer.data["statistics"])]
            assert mock_func_2.call_count == 1
            assert mock_func_2.call_args_list == [
                mock.call(observer.data["statistics"], is_centric=True)
            ]
            assert r == {
                "cc_one_half_plot": "cchalfplot",
                "scaling_tables": "return_tables",
            }
