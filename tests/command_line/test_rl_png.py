from __future__ import annotations

import shutil
import subprocess


def test_rl_png_imported_experiments(dials_data, tmp_path):
    data_dir = dials_data("centroid_test_data", pathlib=True)
    experiments_path = data_dir / "imported_experiments.json"
    strong_pickle = data_dir / "strong.pickle"

    result = subprocess.run(
        [shutil.which("dials.rl_png"), experiments_path, strong_pickle],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    for s in (
        "beam_vector",
        "e3",
        "rotation_axis",
        "solution_1",
        "solution_2",
        "solution_3",
    ):
        assert tmp_path.joinpath(f"rl_{s}.png").is_file()


def test_rl_png_refinement_data(dials_data, tmp_path):
    data_dir = dials_data("refinement_test_data", pathlib=True)
    experiments_path = data_dir / "i04-weak.json"
    indexed_pickle = data_dir / "i04-weak.pickle"

    result = subprocess.run(
        [shutil.which("dials.rl_png"), experiments_path, indexed_pickle],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    for s in ("beam_vector", "e3", "rotation_axis", "a", "b", "c"):
        assert tmp_path.joinpath(f"rl_{s}.png").is_file()
