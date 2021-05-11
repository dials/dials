import os

import procrunner


def test_rl_png_imported_experiments(dials_data, tmpdir):
    data_dir = dials_data("centroid_test_data")
    experiments_path = data_dir / "imported_experiments.json"
    strong_pickle = data_dir / "strong.pickle"

    result = procrunner.run(
        ["dials.rl_png", experiments_path, strong_pickle], working_directory=tmpdir
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
        assert tmpdir.join(f"rl_{s}.png").check()


def test_rl_png_refinement_data(dials_regression, tmpdir):
    data_dir = os.path.join(dials_regression, "refinement_test_data", "i04_weak_data")
    experiments_path = os.path.join(data_dir, "experiments.json")
    indexed_pickle = os.path.join(data_dir, "indexed_strong.pickle")

    result = procrunner.run(
        ["dials.rl_png", experiments_path, indexed_pickle], working_directory=tmpdir
    )
    assert not result.returncode and not result.stderr

    for s in ("beam_vector", "e3", "rotation_axis", "a", "b", "c"):
        assert tmpdir.join(f"rl_{s}.png").check()
