from __future__ import annotations

import json
import logging
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest

from dials.command_line import countrate


def test_generate_histogram_from_shoeboxes():
    shoeboxes = [
        SimpleNamespace(
            data=SimpleNamespace(as_numpy_array=lambda: np.array([[1, 2], [1, -1]]))
        ),
        SimpleNamespace(
            data=SimpleNamespace(as_numpy_array=lambda: np.array([[2, 3]]))
        ),
    ]

    assert countrate.generate_histogram_from_shoeboxes(shoeboxes) == {
        -1: 1,
        1: 2,
        2: 2,
        3: 1,
    }


def test_save_results_to_json(tmp_path):
    results_path = tmp_path / "countrate.json"

    countrate.save_results_to_json({10: 2, 20: 1}, 15302, results_path, 0.25)

    assert json.loads(results_path.read_text()) == {
        "counts": {"10": 2, "20": 1},
        "trusted_range_upper_limit": 15302,
        "recommended_transmission": 0.25,
    }


@pytest.mark.parametrize(
    ("histogram", "percentile", "expected"),
    [
        ({10: 2, 20: 3}, 0, 10),
        ({10: 2, 20: 3}, 0.4, 10),
        ({10: 2, 20: 3}, 0.41, 20),
        ({10: 2, 20: 3}, 1, 20),
    ],
)
def test_get_percentile_intensity(histogram, percentile, expected):
    assert countrate.get_percentile_intensity(histogram, percentile) == expected


@pytest.mark.parametrize("percentile", [-0.01, 1.01])
def test_get_percentile_intensity_rejects_out_of_range_percentile(percentile):
    with pytest.raises(ValueError, match="percentile must be between 0 and 1"):
        countrate.get_percentile_intensity({10: 1}, percentile)


def test_get_percentile_intensity_rejects_empty_histogram():
    with pytest.raises(ValueError, match="empty histogram"):
        countrate.get_percentile_intensity({}, 0.5)


def test_recommend_transmission_scales_transmission():
    assert countrate.recommend_transmission(0.5, 100, 20, 10) == 0.25


def test_recommend_transmission_caps_at_one():
    assert countrate.recommend_transmission(0.8, 100, 10, 50) == 1.0


def _configure_run(
    monkeypatch,
    tmp_path,
    *,
    reflections=None,
    experiments=None,
    has_input=True,
):
    if reflections is None:
        shoeboxes = [
            SimpleNamespace(
                data=SimpleNamespace(as_numpy_array=lambda: np.array([10, 20, 30]))
            )
        ]
        reflections = [{"shoebox": shoeboxes}]
    if experiments is None:
        experiments = [
            SimpleNamespace(
                beam=SimpleNamespace(get_transmission=lambda: 0.5),
                detector=SimpleNamespace(
                    to_dict=lambda: {"panels": [{"trusted_range": [0, 100]}]}
                ),
            )
        ]

    params = SimpleNamespace(
        input=SimpleNamespace(
            experiments=["experiment.expt"] if has_input else [],
            reflections=["reflections.refl"] if has_input else [],
            target_countrate_pct=10.0,
            ref_percentile=90.0,
            trusted_range_warning_threshold=80.0,
        ),
        output=SimpleNamespace(
            log=tmp_path / "countrate.log", results=tmp_path / "countrate.json"
        ),
    )
    options = SimpleNamespace(verbose=0)
    parser = SimpleNamespace(
        parse_args=lambda **kwargs: (params, options),
        diff_phil=SimpleNamespace(as_str=lambda: ""),
        print_help=Mock(),
    )

    monkeypatch.setattr(countrate, "ArgumentParser", lambda **kwargs: parser)
    monkeypatch.setattr(
        countrate,
        "reflections_and_experiments_from_files",
        lambda *args: (reflections, experiments),
    )
    monkeypatch.setattr(countrate.log, "config", lambda **kwargs: None)
    monkeypatch.setattr(countrate, "dials_version", lambda: "test version")
    return parser, params


def test_run_processes_inputs_and_writes_results(monkeypatch, tmp_path, caplog):
    _configure_run(monkeypatch, tmp_path)
    caplog.set_level(logging.INFO, logger="dials.commandline.countrate")

    countrate.run(args=[])

    assert json.loads((tmp_path / "countrate.json").read_text()) == {
        "counts": {"10": 1, "20": 1, "30": 1},
        "trusted_range_upper_limit": 100,
        "recommended_transmission": pytest.approx(1 / 6),
    }
    assert "Found 3 pixels from 1 reflections" in caplog.text
    assert "Recommended max transmission of 16.67%" in caplog.text


def test_run_rejects_missing_input(monkeypatch, tmp_path):
    parser, _ = _configure_run(monkeypatch, tmp_path, has_input=False)

    with pytest.raises(SystemExit, match="Missing required input file"):
        countrate.run(args=[])

    parser.print_help.assert_called_once()


def test_run_requires_exactly_one_experiment_and_reflection_file(monkeypatch, tmp_path):
    parser, _ = _configure_run(monkeypatch, tmp_path, reflections=[], experiments=[])

    with pytest.raises(SystemExit, match="requires exactly one experiment file"):
        countrate.run(args=[])

    parser.print_help.assert_called_once()


def test_analyse_countrate():
    summary = countrate.analyse_countrate(
        histogram={10: 90, 50: 10},
        max_trusted_counts=100,
        transmission=0.5,
        target_countrate_pct=10,
        ref_percentile=90,
        trusted_range_warning_threshold=80,
    )

    assert summary.total_pixels == 100
    assert summary.max_pixel_trusted_range_pct == 50
    assert summary.percentile_trusted_range_pcts == (
        (99.999, 50),
        (99.99, 50),
        (99.9, 50),
        (99.0, 50),
        (90.0, 10),
    )
    assert summary.recommended_transmission == 0.5
    assert summary.ref_percentile_trusted_range_pct == 10
    assert summary.warning_message is None


@pytest.mark.parametrize(
    ("intensity", "expected_warning"),
    [
        (
            50,
            "Some pixels are in the non-linear countrate-correction regime",
        ),
        (
            90,
            "WARNING: The 100th percentile pixel is at 90.00%",
        ),
    ],
)
def test_analyse_countrate_warning(intensity, expected_warning):
    summary = countrate.analyse_countrate(
        histogram={intensity: 100},
        max_trusted_counts=100,
        transmission=0.5,
        target_countrate_pct=10,
        ref_percentile=100,
        trusted_range_warning_threshold=80,
    )

    assert summary.warning_message is not None
    assert expected_warning in summary.warning_message


@pytest.mark.parametrize(
    ("histogram", "trusted_counts", "message"),
    [
        ({}, 100, "empty histogram"),
        ({10: 1}, 0, "max_trusted_counts must be greater than zero"),
    ],
)
def test_analyse_countrate_invalid_input(histogram, trusted_counts, message):
    with pytest.raises(ValueError, match=message):
        countrate.analyse_countrate(
            histogram=histogram,
            max_trusted_counts=trusted_counts,
            transmission=0.5,
            target_countrate_pct=10,
            ref_percentile=90,
            trusted_range_warning_threshold=80,
        )
