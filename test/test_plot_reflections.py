import procrunner


def test_run(dials_data, tmp_path):
    procrunner.run(
        (
            "dials.plot_reflections",
            dials_data("centroid_test_data") / "experiments.json",
            dials_data("centroid_test_data") / "integrated.refl",
            "scan_range=0,5",
        ),
        working_directory=tmp_path,
    ).check_returncode()
    assert (tmp_path / "centroids.png").is_file()
