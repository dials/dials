from __future__ import annotations

import random

import pytest

from dxtbx.serialize import load

from dials.array_family import flex
from dials.util.reciprocal_lattice import Render3d


@pytest.fixture
def multi_sequence_data(dials_data):
    data_dir = dials_data("indexing_test_data", pathlib=True)
    experiments = load.experiment_list(
        data_dir / "multi_sweep-experiments.json",
        check_format=False,
    )
    reflections = flex.reflection_table.from_file(
        data_dir / "multi_sweep-indexed.pickle",
    )
    return {"reflections": reflections, "experiments": experiments}


@pytest.fixture
def centroid_test_data(dials_data):
    experiments = load.experiment_list(
        dials_data("centroid_test_data", pathlib=True) / "experiments.json",
        check_format=False,
    )
    reflections = flex.reflection_table.from_file(
        dials_data("centroid_test_data", pathlib=True) / "integrated.pickle"
    )
    return {"reflections": reflections, "experiments": experiments}


def test_Render3d(mocker, multi_sequence_data):
    experiments = multi_sequence_data["experiments"]
    reflections = multi_sequence_data["reflections"]
    render = Render3d()
    render.viewer = mocker.Mock()
    mocker.spy(render, "set_beam_centre")
    mocker.spy(render, "map_points_to_reciprocal_space")
    mocker.spy(render, "set_points")
    render.load_models(experiments, reflections)
    assert render.set_beam_centre.call_count == 0
    assert render.map_points_to_reciprocal_space.call_count == 1
    assert render.set_points.call_count == 1
    assert render.viewer.set_points.call_count == 1
    assert render.viewer.set_points.call_args[0][0].size() == 1255
    assert render.viewer.set_colors.call_count == 1
    assert render.viewer.set_palette.call_count == 1
    assert render.viewer.set_reciprocal_lattice_vectors.call_count == 1
    assert list(render.viewer.set_palette.call_args[0][0][:2]) == [
        (1.0, 1.0, 1.0),
        (0.9019607843137255, 0.6235294117647059, 0.0),
    ]

    render.settings.beam_centre = (211, 205)
    render.settings.beam_centre = (41, 97)
    render.load_models(experiments, reflections)
    assert render.set_beam_centre.call_count == 1

    for outlier_display, expected_count in (("outliers", 0), ("inliers", 1255)):
        render.settings.outlier_display = outlier_display
        render.load_models(experiments, reflections)
        assert render.viewer.set_points.call_args[0][0].size() == expected_count

    for display, expected_count in (
        ("indexed", 1255),
        ("unindexed", 0),
        ("integrated", 0),
    ):
        render.settings.display = display
        render.load_models(experiments, reflections)
        assert render.viewer.set_points.call_args[0][0].size() == expected_count

    render.settings.display = "all"
    render.settings.experiment_ids = [0, 2, 3]
    render.load_models(experiments, reflections)
    assert render.viewer.set_points.call_args[0][0].size() == 957

    render.settings.black_background = False
    render.load_models(experiments, reflections)
    assert list(render.viewer.set_palette.call_args[0][0][:2]) == [
        (0.0, 0.0, 0.0),
        (0.0980392156862745, 0.3764705882352941, 1.0),
    ]

    render.settings.reverse_phi = True
    mocker.spy(flex.reflection_table, "map_centroids_to_reciprocal_space")
    render.load_models(experiments, reflections)
    assert flex.reflection_table.map_centroids_to_reciprocal_space.call_count == 1


def test_Render3d_integrated(mocker, centroid_test_data):
    experiments = centroid_test_data["experiments"]
    reflections = centroid_test_data["reflections"]

    render = Render3d()
    render.viewer = mocker.Mock()
    render.load_models(experiments, reflections)
    assert render.viewer.set_points.call_args[0][0].size() == 2269

    render.settings.partiality_min = 0.2
    render.settings.partiality_max = 0.9
    render.load_models(experiments, reflections)
    assert render.viewer.set_points.call_args[0][0].size() == 923

    render = Render3d()
    render.viewer = mocker.Mock()
    render.settings.z_min = 2
    render.settings.z_max = 7
    render.load_models(experiments, reflections)
    assert render.viewer.set_points.call_args[0][0].size() == 1885

    random.seed(42)
    reflections["n_signal"] = flex.size_t(
        random.randint(1, 100) for i in range(len(reflections))
    )

    render = Render3d()
    render.viewer = mocker.Mock()
    render.settings.n_min = None
    render.settings.n_max = None
    render.load_models(experiments, reflections)
    assert render.viewer.set_points.call_args[0][0].size() == 2269

    render = Render3d()
    render.viewer = mocker.Mock()
    render.settings.n_min = 20
    render.settings.n_max = 80
    render.load_models(experiments, reflections)
    expected = 1368
    assert render.viewer.set_points.call_args[0][0].size() == expected
