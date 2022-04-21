from __future__ import annotations

import pytest

from cctbx import crystal

from dials.algorithms.indexing.ssx.analysis import (
    combine_results_dicts,
    generate_html_report,
    generate_plots,
    make_cluster_plots,
    make_summary_table,
)


def generate_test_results_dict(n_lattices=1):
    results = {
        0: [
            {
                "Image": "test_image_001.cbf",
                "n_indexed": 0,
                "n_strong": 100,
            }
        ],  # an unindexed image
        1: [
            {
                "Image": "test_image_002.cbf",
                "n_indexed": 50,
                "n_strong": 200,
                "RMSD_X": 1.0,
                "RMSD_Y": 1.1,
                "RMSD_dPsi": 1.2,
            }
        ],  # an indexed image with one lattice
        2: [
            {
                "Image": "test_image_003.cbf",
                "n_indexed": 30,
                "n_strong": 50,
                "RMSD_X": 0.2,
                "RMSD_Y": 0.4,
                "RMSD_dPsi": 0.6,
            }
        ],  # an indexed image with one lattice
    }
    if n_lattices > 1:
        results[2].append(
            {
                "Image": "test_image_003.cbf",
                "n_indexed": 10,
                "n_strong": 50,
                "RMSD_X": 0.3,
                "RMSD_Y": 0.5,
                "RMSD_dPsi": 0.7,
            }
        )  # an indexed image with two lattices
    return results


@pytest.mark.parametrize("n_lattices", [1, 2])
def test_make_summary_table(n_lattices):
    """Test that the summary table has the correct columns"""
    results = generate_test_results_dict(n_lattices)
    table = make_summary_table(results)
    headerline = table.splitlines()[1]
    headers = ["Image", "expt_id", "n_indexed", "RMSD X", "RMSD Y", "RMSD dPsi"]
    last_lattice_line = table.splitlines()[-2]
    assert all(h in headerline for h in headers)
    if n_lattices > 1:
        assert "lattice" in headerline
        expected_last = "| test_image_003.cbf |         2 |         2 | 10/50 (20.0%)  |      0.3 |      0.5 |         0.7 |"
    else:
        assert "lattice" not in headerline
        expected_last = "| test_image_003.cbf |         1 | 30/50 (60.0%)  |      0.2 |      0.4 |         0.6 |"
    assert expected_last == last_lattice_line


def test_combine_multiple_results_summary_dicts():
    s1 = generate_test_results_dict()
    s2 = generate_test_results_dict(n_lattices=2)
    combined = combine_results_dicts([s1, s2])
    assert list(combined.keys()) == [0, 1, 2, 3, 4, 5]
    for k, res in combined.items():
        if k == 5:
            assert len(res) == 2
        else:
            assert len(res) == 1


@pytest.mark.parametrize("n_lattices", [1, 2])
def test_generate_plots(n_lattices):
    plots = generate_plots(generate_test_results_dict(n_lattices))

    # First plot, n_indexed vs image, should have a scatter plot n_strong,
    # one for the first lattice an another for second lattices
    assert len(plots["n_indexed"]["data"]) == n_lattices + 1
    assert len(plots["n_indexed"]["data"][0]["x"]) == 3  # lattice 1 (0 or indexed)
    assert plots["n_indexed"]["data"][0]["y"][0] == 0.0  # first has none indexed
    assert all(plots["n_indexed"]["data"][0]["y"][i] > 0.0 for i in range(1, 3))
    assert plots["n_indexed"]["data"][-1]["y"] == [100, 200, 50]  # n_strong

    if n_lattices == 2:
        assert (
            len(plots["n_indexed"]["data"][1]["x"]) == 1
        )  # only one has second lattice
        assert "lattice 1" in plots["n_indexed"]["data"][0]["name"]
        assert "lattice 2" in plots["n_indexed"]["data"][1]["name"]
        percent_indexed = [0.0, 25.0, 80.0]
    else:
        assert "lattice 1" not in plots["n_indexed"]["data"][0]["name"]
        percent_indexed = [0.0, 25.0, 60.0]

    # Next plot, percentage of strong spots, should just have one trend
    assert plots["percent_indexed"]["data"][0]["y"] == percent_indexed

    # Perecent indexed histogram
    assert sum(plots["percent_indexed_hist"]["data"][0]["y"]) == 3

    # rmsd plots
    assert plots["rmsds"]["data"][0]["y"] == [1.0, 0.2]  # X
    assert plots["rmsds"]["data"][1]["y"] == [1.1, 0.4]  # Y
    assert plots["rmsdz"]["data"][0]["y"] == [1.2, 0.6]  # dPsi
    if n_lattices == 2:
        assert plots["rmsds"]["data"][2]["y"] == [0.3]
        assert plots["rmsds"]["data"][3]["y"] == [0.5]
        assert plots["rmsdz"]["data"][1]["y"] == [0.7]  # dPsi

    # rmsd distribution plots
    assert sum(plots["rmsdxy_hist"]["data"][0]["y"]) == n_lattices + 1
    assert len(plots["rmsdxy_hist"]["data"]) == 2
    assert sum(plots["rmsdz_hist"]["data"][0]["y"]) == n_lattices + 1


def test_generate_html_report(tmp_path):
    plots = generate_plots(generate_test_results_dict())
    fname = "test_report_name.html"
    generate_html_report(plots, tmp_path / fname)
    assert tmp_path.joinpath("test_report_name.html").is_file()


def test_make_cluster_plots():
    from dials.algorithms.clustering.unit_cell import Cluster

    c1 = Cluster(
        [
            crystal.symmetry(
                unit_cell=(10.0, 10.0, 10.0, 90, 90, 90), space_group="P1"
            ),
            crystal.symmetry(
                unit_cell=(10.1, 10.1, 10.1, 90, 90, 90), space_group="P1"
            ),
            crystal.symmetry(
                unit_cell=(10.2, 10.2, 10.2, 90, 90, 90), space_group="P1"
            ),
        ]
    )
    c2 = Cluster(
        [
            crystal.symmetry(
                unit_cell=(11.0, 11.0, 11.0, 90, 90, 90), space_group="P1"
            ),
            crystal.symmetry(
                unit_cell=(11.1, 11.1, 11.1, 90, 90, 90), space_group="P1"
            ),
            crystal.symmetry(
                unit_cell=(11.2, 11.2, 11.2, 90, 90, 90), space_group="P1"
            ),
            crystal.symmetry(
                unit_cell=(11.3, 11.3, 11.3, 90, 90, 90), space_group="P1"
            ),
        ]
    )
    clusters = [c1, c2]
    plots = make_cluster_plots(clusters)
    assert "uc_scatter_0" in plots
    assert "uc_scatter_1" in plots
    assert "uc_hist_0" in plots
    assert "uc_hist_1" in plots
    print(plots)
    assert len(plots["uc_hist_0"]["data"]) == 3
    assert len(plots["uc_hist_0"]["data"][0]["x"]) == 3
    assert len(plots["uc_hist_1"]["data"][0]["x"]) == 4
    assert len(plots["uc_scatter_0"]["data"]) == 3
    assert len(plots["uc_scatter_0"]["data"][0]["x"]) == 3
    assert len(plots["uc_scatter_1"]["data"][0]["x"]) == 4
