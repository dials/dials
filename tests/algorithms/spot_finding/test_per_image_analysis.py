from __future__ import annotations

import pytest

from dxtbx.serialize import load

from dials.algorithms.spot_finding import per_image_analysis
from dials.array_family import flex


@pytest.fixture
def centroid_test_data(dials_data):
    reflections = flex.reflection_table.from_file(
        dials_data("centroid_test_data", pathlib=True) / "strong.pickle"
    )
    experiments = load.experiment_list(
        dials_data("centroid_test_data", pathlib=True) / "imported_experiments.json"
    )
    reflections.centroid_px_to_mm(experiments)
    reflections.map_centroids_to_reciprocal_space(experiments)
    return experiments, reflections


def test_estimate_resolution_limit_distl_method1(centroid_test_data, tmp_path):
    experiments, reflections = centroid_test_data
    plot_file = tmp_path / "distl1.png"
    d_min, noisiness = per_image_analysis.estimate_resolution_limit_distl_method1(
        reflections=reflections, plot_filename=plot_file
    )
    assert d_min == pytest.approx(1.8771983880778702)
    assert noisiness == pytest.approx(0.021021021021021023)
    assert plot_file.is_file()


def test_estimate_resolution_limit_distl_method2(centroid_test_data, tmp_path):
    experiments, reflections = centroid_test_data
    plot_file = tmp_path / "distl2.png"
    d_min, noisiness = per_image_analysis.estimate_resolution_limit_distl_method2(
        reflections=reflections, plot_filename=plot_file
    )
    assert d_min == pytest.approx(1.6293601446185495)
    assert noisiness == pytest.approx(0.0858974358974359)
    assert plot_file.is_file()


def test_estimate_resolution(centroid_test_data, tmp_path):
    experiments, reflections = centroid_test_data
    ice_sel = per_image_analysis.ice_rings_selection(reflections, width=0.004)
    assert ice_sel.count(True) == 76
    plot_file = tmp_path / "i_over_sigi_vs_resolution.png"
    d_min = per_image_analysis.estimate_resolution_limit(
        reflections=reflections, plot_filename=plot_file
    )
    assert d_min == pytest.approx(1.446715534174674)
    assert plot_file.is_file()
    d_min = per_image_analysis.estimate_resolution_limit(
        reflections=reflections, ice_sel=ice_sel
    )
    assert d_min == pytest.approx(1.446715534174674)


def test_stats_for_reflection_table(centroid_test_data):
    _, reflections = centroid_test_data
    stats = per_image_analysis.stats_for_reflection_table(reflections)
    result = stats._asdict()
    expected = {
        "d_min_distl_method_1": 1.8771983880778702,
        "d_min_distl_method_2": 1.6293601446185495,
        "estimated_d_min": 1.446715534174674,
        "n_spots_4A": 76,
        "n_spots_no_ice": 578,
        "n_spots_total": 654,
        "noisiness_method_1": 0.021021021021021023,
        "noisiness_method_2": 0.0858974358974359,
        "total_intensity": 919847.0,
    }
    for k, v in expected.items():
        assert result[k] == pytest.approx(v)


def test_stats_for_reflection_table_no_resolution_analysis_no_ice_filtering(
    centroid_test_data,
):
    _, reflections = centroid_test_data
    stats = per_image_analysis.stats_for_reflection_table(
        reflections, resolution_analysis=False, filter_ice=False
    )
    result = stats._asdict()
    expected = {
        "d_min_distl_method_1": -1.0,
        "d_min_distl_method_2": -1.0,
        "estimated_d_min": -1.0,
        "n_spots_4A": 76,
        "n_spots_no_ice": 654,
        "n_spots_total": 654,
        "noisiness_method_1": -1.0,
        "noisiness_method_2": -1.0,
        "total_intensity": 961472.0,
    }
    for k, v in expected.items():
        assert result[k] == pytest.approx(v)


def test_stats_per_image(centroid_test_data):
    experiments, reflections = centroid_test_data
    stats = per_image_analysis.stats_per_image(experiments[0], reflections)
    result = stats._asdict()
    for v in result.values():
        assert len(v) == len(experiments[0].scan)
    assert stats.n_spots_total == [90, 100, 67, 49, 54, 62, 68, 83, 81]
    t = stats.as_table()
    assert len(t) == len(experiments[0].scan) + 1
    assert t[0] == [
        "image",
        "#spots",
        "#spots_no_ice",
        "total_intensity",
        "d_min",
        "d_min (distl method 1)",
        "d_min (distl method 2)",
    ]
    assert t[1] == ["1", "90", "77", "28214", "1.56", "2.08 (0.09)", "1.59 (0.27)"]
    # Test n_rows option
    t = stats.as_table(n_rows=3)
    assert len(t) == 4
    # Test perm option
    perm = flex.random_permutation(len(experiments[0].scan))
    t = stats.as_table(perm=perm)
    assert [tt[0] for tt in t[1:]] == [str(i + 1) for i in perm]


def test_stats_table_no_resolution_analysis(centroid_test_data):
    experiments, reflections = centroid_test_data
    stats = per_image_analysis.stats_per_image(
        experiments[0], reflections, resolution_analysis=False
    )
    t = stats.as_table()
    assert t[0] == ["image", "#spots", "#spots_no_ice", "total_intensity"]


def test_plot_stats(centroid_test_data, tmp_path):
    experiments, reflections = centroid_test_data
    stats = per_image_analysis.stats_per_image(
        experiments[0], reflections, resolution_analysis=False
    )
    image_file = tmp_path / "pia.png"
    per_image_analysis.plot_stats(stats, filename=image_file)
    assert image_file.is_file()
