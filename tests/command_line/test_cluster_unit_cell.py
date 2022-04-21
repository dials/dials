from __future__ import annotations

import glob
import os
import pathlib

import procrunner
import pytest

from cctbx import crystal
from dxtbx.model.experiment_list import ExperimentList, ExperimentListFactory
from dxtbx.serialize import load

from dials.command_line import cluster_unit_cell


def test_dials_cluster_unit_cell_command_line(dials_regression, tmp_path):
    pytest.importorskip("scipy")
    pytest.importorskip("xfel")

    data_dir = (
        pathlib.Path(dials_regression) / "refinement_test_data" / "multi_narrow_wedges"
    )
    experiments = list(data_dir.glob("data/sweep_*/experiments.json"))

    result = procrunner.run(
        command=["dials.cluster_unit_cell", "plot.show=False"] + experiments,
        print_stdout=False,
        working_directory=tmp_path,
    )
    assert not result.returncode
    assert tmp_path.joinpath("cluster_unit_cell.png").is_file()


def test_dials_cluster_unit_cell_command_line_output_files(dials_regression, tmp_path):
    pytest.importorskip("scipy")
    pytest.importorskip("xfel")

    data_dir = (
        pathlib.Path(dials_regression) / "refinement_test_data" / "multi_narrow_wedges"
    )
    experiments = list(data_dir.glob("data/sweep_*/experiments.json"))
    reflections = list(data_dir.glob("data/sweep_*/reflections.pickle"))

    # Combine experiments. Write PHIL file to avoid "command line is too long" error on Windows
    with open(tmp_path / "input.phil", "w") as f:
        f.writelines((f"input.reflections={i}" + os.linesep for i in reflections))
        f.writelines((f"input.experiments={i}" + os.linesep for i in experiments))
    result = procrunner.run(
        command=["dials.combine_experiments", "input.phil"],
        working_directory=tmp_path,
    )

    assert not result.returncode
    assert (tmp_path / "combined.refl").is_file()
    assert (tmp_path / "combined.expt").is_file()

    result = procrunner.run(
        command=[
            "dials.cluster_unit_cell",
            "plot.show=False",
            tmp_path / "combined.refl",
            tmp_path / "combined.expt",
            "output.clusters=True",
            "threshold=40",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode
    assert (tmp_path / "cluster_unit_cell.png").is_file()
    assert (tmp_path / "cluster_0.refl").is_file()
    assert (tmp_path / "cluster_0.expt").is_file()
    expts = load.experiment_list(tmp_path / "cluster_0.expt", check_format=False)
    assert len(expts) == 101
    assert (tmp_path / "cluster_1.refl").is_file()
    assert (tmp_path / "cluster_1.expt").is_file()
    expts = load.experiment_list(tmp_path / "cluster_1.expt", check_format=False)
    assert len(expts) == 1
    assert (tmp_path / "cluster_2.refl").is_file()
    assert (tmp_path / "cluster_2.expt").is_file()
    expts = load.experiment_list(tmp_path / "cluster_2.expt", check_format=False)
    assert len(expts) == 1

    result = procrunner.run(
        command=[
            "dials.split_experiments",
            tmp_path / "combined.refl",
            tmp_path / "combined.expt",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode
    experiments = list(tmp_path.glob("split_*.expt"))
    reflections = list(tmp_path.glob("split_*.refl"))

    # Write PHIL file to avoid "command line is too long" error on Windows
    with open(tmp_path / "input.phil", "w") as f:
        f.writelines((f"input.reflections={i}" + os.linesep for i in reflections))
        f.writelines((f"input.experiments={i}" + os.linesep for i in experiments))
    result = procrunner.run(
        command=[
            "dials.cluster_unit_cell",
            "output.clusters=True",
            "threshold=40",
            "plot.show=False",
            "input.phil",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode


def test_cluster_unit_cell_api(dials_regression):
    pytest.importorskip("scipy")
    pytest.importorskip("xfel")

    data_dir = os.path.join(
        dials_regression, "refinement_test_data", "multi_narrow_wedges"
    )
    experiments = ExperimentList(
        [
            ExperimentListFactory.from_json_file(expt, check_format=False)[0]
            for expt in glob.glob(
                os.path.join(data_dir, "data/sweep_*/experiments.json")
            )
        ]
    )
    crystal_symmetries = [
        crystal.symmetry(
            unit_cell=expt.crystal.get_unit_cell(),
            space_group=expt.crystal.get_space_group(),
        )
        for expt in experiments
    ]

    params = cluster_unit_cell.phil_scope.extract()
    params.plot.show = False
    params.plot.name = None
    clusters = cluster_unit_cell.do_cluster_analysis(crystal_symmetries, params)
    assert len(clusters) == 1
    cluster = clusters[0]
    assert len(cluster) == 40
    assert cluster.median_cell == pytest.approx(
        [
            90.9430182020995,
            90.9430182020995,
            90.9430182020995,
            109.47122063449069,
            109.47122063449069,
            109.47122063449069,
        ],
        abs=1e-6,
    )
    assert cluster.cell_std == pytest.approx(
        [0.09509739126548639, 0.09509739126548526, 0.0950973912654865, 0, 0, 0],
        abs=1e-6,
    )
