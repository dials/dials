from __future__ import annotations

import pickle
import shutil
import subprocess

import pytest

from dxtbx.model import ExperimentList
from scitbx.array_family import flex


@pytest.mark.skip(reason="Apparently causes SEGV on some platforms")
def test_model_background(dials_data, tmp_path):
    # Use a data set from a P2M for speed (small detector).
    data_dir = dials_data("l_cysteine_dials_output", pathlib=True)

    result = subprocess.run(
        [shutil.which("dials.model_background"), data_dir / "11_integrated.expt"],
        cwd=tmp_path,
        capture_output=True,
    )
    print(result.stderr.decode())
    result.check_returncode()
    assert not result.stderr
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
    assert data.all() == (1679, 1475)
    min_max_mean = flex.min_max_mean_double(data.as_1d())
    assert min_max_mean.max == pytest.approx(0.1800928143817288)
    assert min_max_mean.min == 0.0
    assert min_max_mean.mean == pytest.approx(0.020816338853321865)

    # Test integration using this background model. It turns out that 11_integrated.{expt,refl}
    # can't be passed to dials.integrate, so we will make our own input from the
    # indexed.{expt,refl}. This has 4 experiments and 1700 images, and we'll take just the first
    # experiment and 10 images for this test.
    subprocess.run(
        [
            shutil.which("dials.split_experiments"),
            data_dir / "indexed.expt",
            data_dir / "indexed.refl",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    expts = ExperimentList.from_file(tmp_path / "split_0.expt")
    expts[0].scan.set_image_range((1, 10))
    expts.as_file(tmp_path / "modified.expt")

    # Test integration using the background model, with robust.algorithm=(True|False)
    result = subprocess.run(
        [
            shutil.which("dials.integrate"),
            tmp_path / "modified.expt",
            tmp_path / "split_0.refl",
            "background.algorithm=gmodel",
            "gmodel.robust.algorithm=False",
            "gmodel.model=background.pickle",
        ],
        cwd=tmp_path,
        capture_output=True,
    )

    assert not result.returncode and not result.stderr

    result = subprocess.run(
        [
            shutil.which("dials.integrate"),
            tmp_path / "modified.expt",
            tmp_path / "split_0.refl",
            "background.algorithm=gmodel",
            "gmodel.robust.algorithm=True",
            "gmodel.model=background.pickle",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
