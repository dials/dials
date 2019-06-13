from __future__ import absolute_import, division, print_function

import os
import pytest

from dxtbx.serialize import load
from dials.array_family import flex
from dials.util.reciprocal_lattice import Render3d


def test_Render3d(mocker, dials_regression):
    experiments = load.experiment_list(
        os.path.join(
            dials_regression, "indexing_test_data", "multi_sweep", "experiments.json"
        ),
        check_format=False,
    )
    reflections = flex.reflection_table.from_file(
        os.path.join(
            dials_regression, "indexing_test_data", "multi_sweep", "indexed.pickle"
        )
    )
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

    for (outlier_display, expected_count) in (("outliers", 0), ("inliers", 1254)):
        render.settings.outlier_display = outlier_display
        render.load_models(experiments, reflections)
        assert render.viewer.set_points.call_args[0][0].size() == expected_count

    for (display, expected_count) in (
        ("indexed", 1254),
        ("unindexed", 0),
        ("integrated", 0),
    ):
        render.settings.display = display
        render.load_models(experiments, reflections)
        assert render.viewer.set_points.call_args[0][0].size() == expected_count

    render.settings.display = "all"
    render.settings.experiment_ids = [0, 2, 3]
    render.load_models(experiments, reflections)
    assert render.viewer.set_points.call_args[0][0].size() == 956

    render.settings.black_background = False
    render.load_models(experiments, reflections)
    assert list(render.viewer.set_palette.call_args[0][0][:2]) == [
        (0.0, 0.0, 0.0),
        (0.0980392156862745, 0.3764705882352941, 1.0),
    ]

    render.settings.reverse_phi = True
    mocker.spy(flex.reflection_table, "map_centroids_to_reciprocal_space")
    render.load_models(experiments, reflections)
    assert flex.reflection_table.map_centroids_to_reciprocal_space.call_count == len(
        experiments
    )
    assert flex.reflection_table.map_centroids_to_reciprocal_space.call_args[0][
        3
    ].get_rotation_axis() == pytest.approx((-1.0, 0.0, 0.0))
