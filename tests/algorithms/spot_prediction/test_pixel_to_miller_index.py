from __future__ import annotations

import pytest


def test(dials_data):
    from dxtbx.model.experiment_list import ExperimentListFactory

    from dials.algorithms.spot_prediction import PixelToMillerIndex
    from dials.array_family import flex

    filename = dials_data("centroid_test_data", pathlib=True) / "experiments.json"

    experiments = ExperimentListFactory.from_json_file(filename)

    transform = PixelToMillerIndex(
        experiments[0].beam,
        experiments[0].detector,
        experiments[0].goniometer,
        experiments[0].scan,
        experiments[0].crystal,
    )

    reflections = flex.reflection_table.from_predictions(experiments[0])

    for r in reflections.rows():
        panel = r["panel"]
        x, y, z = r["xyzcal.px"]
        h0 = r["miller_index"]
        h1 = transform.h(panel, x, y, z)
        assert h0 == pytest.approx(h1, abs=1e-7)
