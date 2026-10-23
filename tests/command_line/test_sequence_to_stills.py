from __future__ import annotations

import shutil
import subprocess

from dxtbx.model.experiment_list import ExperimentListFactory


def test_sequence_to_stills(dials_data, tmp_path):
    data_dir = dials_data("insulin_processed")
    input_experiments = data_dir / "integrated.expt"
    input_reflections = data_dir / "refined.refl"
    result = subprocess.run(
        [
            shutil.which("dials.sequence_to_stills"),
            input_experiments,
            input_reflections,
            "domain_size_ang=500",
            "half_mosaicity_deg=0.1",
            "max_scan_points=10",
        ],
        cwd=tmp_path,
    )
    assert not result.returncode and not result.stderr

    assert (tmp_path / "stills.expt").is_file()
    assert (tmp_path / "stills.refl").is_file()

    experiments = ExperimentListFactory.from_json_file(
        tmp_path / "stills.expt", check_format=False
    )
    assert len(experiments) == 10
    assert len(experiments.identifiers()) == 10


def test_data_with_static_model(dials_data, tmp_path):
    # test for regression of https://github.com/dials/dials/issues/2516
    data_dir = dials_data("insulin_processed")
    input_experiments = data_dir / "indexed.expt"
    input_reflections = data_dir / "indexed.refl"
    result = subprocess.run(
        [
            shutil.which("dials.sequence_to_stills"),
            input_experiments,
            input_reflections,
        ],
        cwd=tmp_path,
    )
    assert not result.returncode and not result.stderr

    assert (tmp_path / "stills.expt").is_file()
    assert (tmp_path / "stills.refl").is_file()

    experiments = ExperimentListFactory.from_json_file(
        tmp_path / "stills.expt", check_format=False
    )
    assert len(experiments) == 45


def test_sliced_sequence(dials_data, tmp_path):
    # test for regression of https://github.com/dials/dials/issues/2519
    data_dir = dials_data("insulin_processed")
    input_experiments = data_dir / "indexed.expt"
    input_reflections = data_dir / "indexed.refl"

    result = subprocess.run(
        [
            shutil.which("dials.slice_sequence"),
            input_experiments,
            input_reflections,
            "image_range=5,45",
        ],
        cwd=tmp_path,
    )
    assert not result.returncode and not result.stderr

    result = subprocess.run(
        [
            shutil.which("dials.sequence_to_stills"),
            "indexed_5_45.expt",
            "indexed_5_45.refl",
        ],
        cwd=tmp_path,
    )
    assert not result.returncode and not result.stderr

    assert (tmp_path / "stills.expt").is_file()
    assert (tmp_path / "stills.refl").is_file()

    experiments = ExperimentListFactory.from_json_file(
        tmp_path / "stills.expt", check_format=False
    )
    assert len(experiments) == 41


def test_shoeboxes_without_background(dials_data):
    # Shoeboxes from spot finding have no background allocated, meaning zero.
    # They must give the same stills as the same shoeboxes with a background of
    # zeros, as they have after a round trip through a reflection file.
    from dials.array_family import flex
    from dials.command_line.sequence_to_stills import phil_scope, sequence_to_stills
    from dials.model.data import Shoebox

    data_dir = dials_data("insulin_processed")
    experiments = ExperimentListFactory.from_json_file(
        data_dir / "integrated.expt", check_format=False
    )
    reflections = flex.reflection_table.from_file(data_dir / "refined.refl")
    reflections = reflections.select(reflections["bbox"].parts()[4] < 5)
    assert all(s.background.all_eq(0) for s in reflections["shoebox"])
    params = phil_scope.extract()
    params.output.domain_size_ang = 500
    params.output.half_mosaicity_deg = 0.1
    params.max_scan_points = 5

    without = reflections.copy()
    shoeboxes = flex.shoebox()
    for s in reflections["shoebox"]:
        new_sb = Shoebox(s.panel, s.bbox)
        new_sb.data = s.data
        new_sb.mask = s.mask
        shoeboxes.append(new_sb)
    without["shoebox"] = shoeboxes

    _, expected = sequence_to_stills(experiments, [reflections], params)
    _, result = sequence_to_stills(experiments, [without], params)
    assert len(result) == len(expected) > 0
    for key in ("intensity.sum.value", "intensity.sum.variance", "xyzobs.px.value"):
        assert list(result[key]) == list(expected[key])
    assert not any(result["shoebox"].is_background_allocated())
    assert all(result["shoebox"].is_consistent())
