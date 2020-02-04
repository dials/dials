from __future__ import absolute_import, division, print_function

import pytest

from dxtbx.serialize import load
from dials.algorithms.spot_finding import per_image_analysis
from dials.array_family import flex


def test_estimate_resolution_limit_distl_method1(dials_data, tmpdir):
    reflections = flex.reflection_table.from_file(
        dials_data("centroid_test_data").join("strong.pickle").strpath
    )
    experiments = load.experiment_list(
        dials_data("centroid_test_data").join("imported_experiments.json")
    )
    reflections.centroid_px_to_mm(experiments)
    reflections.map_centroids_to_reciprocal_space(experiments)
    plot_file = tmpdir.join("distl1.png")
    d_min, noisiness = per_image_analysis.estimate_resolution_limit_distl_method1(
        reflections=reflections, plot_filename=plot_file.strpath
    )
    assert d_min == pytest.approx(1.8771983880778702)
    assert noisiness == pytest.approx(0.021021021021021023)
    assert plot_file.check()


def test_estimate_resolution_limit_distl_method2(dials_data, tmpdir):
    reflections = flex.reflection_table.from_file(
        dials_data("centroid_test_data").join("strong.pickle").strpath
    )
    experiments = load.experiment_list(
        dials_data("centroid_test_data").join("imported_experiments.json")
    )
    reflections.centroid_px_to_mm(experiments)
    reflections.map_centroids_to_reciprocal_space(experiments)
    plot_file = tmpdir.join("distl2.png")
    d_min, noisiness = per_image_analysis.estimate_resolution_limit_distl_method2(
        reflections=reflections, plot_filename=plot_file.strpath
    )
    assert d_min == pytest.approx(1.6293601446185495)
    assert noisiness == pytest.approx(0.0858974358974359)
    assert plot_file.check()


def test_estimate_resolution(dials_data, tmpdir):
    reflections = flex.reflection_table.from_file(
        dials_data("centroid_test_data").join("strong.pickle").strpath
    )
    experiments = load.experiment_list(
        dials_data("centroid_test_data").join("imported_experiments.json")
    )
    reflections.centroid_px_to_mm(experiments)
    reflections.map_centroids_to_reciprocal_space(experiments)
    ice_sel = per_image_analysis.ice_rings_selection(reflections, width=0.004)
    assert ice_sel.count(True) == 76
    plot_file = tmpdir.join("i_over_sigi_vs_resolution.png")
    d_min = per_image_analysis.estimate_resolution_limit(
        reflections=reflections, plot_filename=plot_file.strpath
    )
    assert d_min == pytest.approx(1.446715534174674)
    assert plot_file.check()
    d_min = per_image_analysis.estimate_resolution_limit(
        reflections=reflections, ice_sel=ice_sel
    )
    assert d_min == pytest.approx(1.446715534174674)
