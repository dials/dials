from __future__ import annotations

import pickle
import shutil
import subprocess
import sys

import pytest

from dxtbx.model import ExperimentList
from scitbx.array_family import flex


def test_model_background(dials_data, tmp_path):
    # Use a data set from a P2M for speed (small detector).
    data_dir = dials_data("l_cysteine_dials_output")

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


# Run dials.model_background with the "spawn" start method, as used on Windows,
# which requires everything sent to the worker processes to be picklable
SPAWN_SCRIPT = """
import multiprocessing
import sys

import libtbx.easy_mp  # noqa: F401 - may set the start method on import

if __name__ == "__main__":
    multiprocessing.set_start_method("spawn", force=True)
    from dials.command_line.model_background import run

    run(sys.argv[1:])
"""


def test_model_background_nproc(dials_data, tmp_path):
    data_dir = dials_data("l_cysteine_dials_output")
    (tmp_path / "spawn.py").write_text(SPAWN_SCRIPT)

    models = []
    for name, command in (
        ("nproc_1", [shutil.which("dials.model_background"), "nproc=1"]),
        ("nproc_2", [shutil.which("dials.model_background"), "nproc=2"]),
        ("spawn", [sys.executable, tmp_path / "spawn.py", "nproc=2"]),
    ):
        working_dir = tmp_path / name
        working_dir.mkdir()
        result = subprocess.run(
            command + [data_dir / "11_integrated.expt"],
            cwd=working_dir,
            capture_output=True,
        )
        result.check_returncode()
        assert not result.stderr
        if name != "nproc_1":
            assert b"Using multiprocessing with 2 parallel job(s)" in result.stdout
        with open(working_dir / "background.pickle", "rb") as f:
            models.append(pickle.load(f).data(0))

    for model in models[1:]:
        assert model.all() == models[0].all()
        assert list(model) == list(models[0])


def test_model_background_nproc_releases_results(dials_data, tmp_path, monkeypatch):
    # The parallel map keeps every job result until all have finished, so each
    # must be released once it has been accumulated to avoid running out of memory
    from dials.algorithms.integration import image_integrator
    from dials.command_line import model_background

    original_parallel_map = image_integrator.multi_node_parallel_map
    returned = []

    def parallel_map(*args, **kwargs):
        result = original_parallel_map(*args, **kwargs)
        returned.extend(result)
        return result

    monkeypatch.setattr(image_integrator, "multi_node_parallel_map", parallel_map)
    monkeypatch.chdir(tmp_path)

    model_background.run(
        [str(dials_data("l_cysteine_dials_output") / "11_integrated.expt"), "nproc=2"]
    )
    assert len(returned) > 1
    assert all(result == [] for result in returned)
