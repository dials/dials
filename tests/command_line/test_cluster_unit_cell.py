from __future__ import annotations

import glob
import os
import shutil
import subprocess

import pytest

from cctbx import crystal
from dxtbx.model.experiment_list import ExperimentList, ExperimentListFactory
from dxtbx.serialize import load

from dials.command_line import cluster_unit_cell


def test_dials_cluster_unit_cell_command_line(dials_data, tmp_path):
    pytest.importorskip("scipy")

    data_dir = dials_data("polyhedra_narrow_wedges", pathlib=True)
    experiments = sorted(data_dir.glob("sweep_*_experiments.json"))

    result = subprocess.run(
        [shutil.which("dials.cluster_unit_cell"), "plot.show=False"] + experiments,
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode
    assert tmp_path.joinpath("cluster_unit_cell.png").is_file()


def test_dials_cluster_unit_cell_command_line_output_files(dials_data, tmp_path):
    pytest.importorskip("scipy")

    data_dir = dials_data("polyhedra_narrow_wedges", pathlib=True)
    experiments = sorted(data_dir.glob("sweep_*_experiments.json"))
    reflections = sorted(data_dir.glob("sweep_*_reflections.pickle"))

    # Combine experiments. Write PHIL file to avoid "command line is too long" error on Windows
    with open(tmp_path / "input.phil", "w") as f:
        f.writelines(f"input.reflections={i}" + os.linesep for i in reflections)
        f.writelines(f"input.experiments={i}" + os.linesep for i in experiments)
    result = subprocess.run(
        [shutil.which("dials.combine_experiments"), "input.phil"],
        cwd=tmp_path,
        capture_output=True,
    )

    assert not result.returncode
    assert (tmp_path / "combined.refl").is_file()
    assert (tmp_path / "combined.expt").is_file()

    result = subprocess.run(
        [
            shutil.which("dials.cluster_unit_cell"),
            "plot.show=False",
            tmp_path / "combined.refl",
            tmp_path / "combined.expt",
            "output.clusters=True",
            "threshold=40",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode
    assert (tmp_path / "cluster_unit_cell.png").is_file()
    assert (tmp_path / "cluster_1.refl").is_file()
    assert (tmp_path / "cluster_1.expt").is_file()
    expts = load.experiment_list(tmp_path / "cluster_1.expt", check_format=False)
    assert len(expts) == 101
    assert (tmp_path / "cluster_2.refl").is_file()
    assert (tmp_path / "cluster_2.expt").is_file()
    expts = load.experiment_list(tmp_path / "cluster_2.expt", check_format=False)
    assert len(expts) == 1
    assert (tmp_path / "cluster_3.refl").is_file()
    assert (tmp_path / "cluster_3.expt").is_file()
    expts = load.experiment_list(tmp_path / "cluster_3.expt", check_format=False)
    assert len(expts) == 1

    result = subprocess.run(
        [
            shutil.which("dials.split_experiments"),
            tmp_path / "combined.refl",
            tmp_path / "combined.expt",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode
    experiments = list(tmp_path.glob("split_*.expt"))
    reflections = list(tmp_path.glob("split_*.refl"))

    # Write PHIL file to avoid "command line is too long" error on Windows
    with open(tmp_path / "input.phil", "w") as f:
        f.writelines(f"input.reflections={i}" + os.linesep for i in reflections)
        f.writelines(f"input.experiments={i}" + os.linesep for i in experiments)
    result = subprocess.run(
        [
            shutil.which("dials.cluster_unit_cell"),
            "output.clusters=True",
            "threshold=40",
            "plot.show=False",
            "input.phil",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode


def test_cluster_unit_cell_api(dials_data):
    pytest.importorskip("scipy")

    data_dir = dials_data("polyhedra_narrow_wedges", pathlib=True)
    experiments = ExperimentList(
        [
            ExperimentListFactory.from_json_file(expt, check_format=False)[0]
            for expt in glob.glob(os.path.join(data_dir, "sweep_*_experiments.json"))
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
