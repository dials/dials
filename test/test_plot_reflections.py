import subprocess


def test_run(dials_data, tmp_path):
    subprocess.run(
        (
            "dials.plot_reflections",
            dials_data("centroid_test_data") / "experiments.json",
            dials_data("centroid_test_data") / "integrated.refl",
            "scan_range=0,5",
        ),
        check=True,
        cwd=tmp_path,
    )
    assert (tmp_path / "centroids.png").is_file()
