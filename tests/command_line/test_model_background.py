from __future__ import annotations

import pickle

import procrunner
import pytest

from scitbx.array_family import flex


def test_model_background(dials_data, tmp_path):
    centroid = dials_data("centroid_test_data", pathlib=True)
    expts = centroid / "experiments.json"

    result = procrunner.run(
        ["dials.model_background", expts],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    for filename in (
        "background.pickle",
        "mean_0.png",
        "variance_0.png",
        "dispersion_0.png",
        "mask_0.png",
        "min_0.png",
        "max_0.png",
        "model_0.png",
    ):
        assert (tmp_path / filename).is_file()

    with open((tmp_path / "background.pickle"), "rb") as f:
        background = pickle.load(f)

    panel = 0
    data = background.data(panel)
    assert data.all() == (2527, 2463)
    min_max_mean = flex.min_max_mean_double(data.as_1d())
    assert min_max_mean.max == pytest.approx(5.9114028830604095)
    assert min_max_mean.min == 0.0
    assert min_max_mean.mean == pytest.approx(0.5013730161480899)

    # Test integration using the background model, with robust.algorithm=(True|False)
    refls = centroid / "indexed.refl"

    result = procrunner.run(
        [
            "dials.integrate",
            expts,
            refls,
            "background.algorithm=gmodel",
            "gmodel.robust.algorithm=False",
            "gmodel.model=background.pickle",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr

    result = procrunner.run(
        [
            "dials.integrate",
            expts,
            refls,
            "background.algorithm=gmodel",
            "gmodel.robust.algorithm=True",
            "gmodel.model=background.pickle",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
