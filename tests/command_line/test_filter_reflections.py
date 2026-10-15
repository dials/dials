from __future__ import annotations

import shutil
import subprocess

import pytest

from dxtbx.util import ersatz_uuid4

from dials.array_family import flex


@pytest.fixture(scope="session")
def reflections(tmp_path_factory):
    # Make a dummy reflection table for the test setting some values and flags
    rt = flex.reflection_table.empty_standard(6)
    rt["iobs"] = flex.size_t_range(len(rt))
    rt["panel"] = flex.size_t_range(len(rt))
    rt["id"] = flex.int([0] * 5 + [1])
    rt["d"] = flex.double([50, 40, 3.0, 2.5, 2.0, 1.0])
    mask1 = flex.bool([True] * 3 + [False] * 3)
    mask2 = flex.bool([True, False] * 3)
    rt.set_flags(mask1, rt.flags.integrated)
    rt.set_flags(mask2, rt.flags.reference_spot)
    tmp_path = tmp_path_factory.mktemp("filter_reflections")
    rt_name = tmp_path / "test_refs.refl"
    for id in set(rt["id"]):
        rt.experiment_identifiers()[id] = ersatz_uuid4()
    rt.as_file(rt_name)
    return rt_name


def test_filter_reflections_flag_expression(reflections, tmp_path):
    result = subprocess.run(
        [
            shutil.which("dials.filter_reflections"),
            reflections,
            "flag_expression='integrated & ~reference_spot'",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    ref = flex.reflection_table.from_file(tmp_path / "filtered.refl")
    # The test selects only the 2nd reflection
    assert len(ref) == 1
    assert list(ref["iobs"]) == [1]


def test_filter_reflections_by_experiment_id(reflections, tmp_path):
    result = subprocess.run(
        [shutil.which("dials.filter_reflections"), reflections, "id=0"],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    ref = flex.reflection_table.from_file(tmp_path / "filtered.refl")
    # The test selects only the first five reflections
    assert len(ref) == 5
    assert list(ref["iobs"]) == [0, 1, 2, 3, 4]


def test_filter_reflections_by_panel(reflections, tmp_path):
    result = subprocess.run(
        [shutil.which("dials.filter_reflections"), reflections, "panel=5"],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    ref = flex.reflection_table.from_file(tmp_path / "filtered.refl")
    # The test selects only the last reflection
    assert len(ref) == 1
    assert list(ref["iobs"]) == [5]


def test_filter_reflections_by_resolution(reflections, tmp_path):
    result = subprocess.run(
        [
            shutil.which("dials.filter_reflections"),
            reflections,
            "d_max=3.0",
            "d_min=2.0",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    ref = flex.reflection_table.from_file(tmp_path / "filtered.refl")
    # The test selects only the 3rd, 4th and 5th reflections
    assert len(ref) == 3
    assert list(ref["iobs"]) == [2, 3, 4]


def test_filter_reflections_printing_analysis(reflections, tmp_path):
    result = subprocess.run(
        [shutil.which("dials.filter_reflections"), reflections],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr


@pytest.mark.parametrize(
    "experiment,reflections,args,expected",
    [
        (
            "imported_experiments.json",
            "strong.pickle",
            ["ice_rings.filter=True"],
            615,
        ),
        (
            "imported_experiments.json",
            "strong.pickle",
            ["d_min=3.1", "d_max=20"],
            153,
        ),
        (
            "experiments.json",
            "integrated.refl",
            ["d_min=2", "d_max=20"],
            371,
        ),
    ],
)
def test_filter_reflections(
    experiment, reflections, args, expected, dials_data, tmp_path
):
    dataset = dials_data("centroid_test_data")
    result = subprocess.run(
        [
            shutil.which("dials.filter_reflections"),
            dataset / experiment,
            dataset / reflections,
        ]
        + args,
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    filtered_refl = flex.reflection_table.from_file(tmp_path / "filtered.refl")
    assert len(filtered_refl) == expected


def test_filter_reflections_by_reflection_index(reflections, tmp_path):
    result = subprocess.run(
        [
            shutil.which("dials.filter_reflections"),
            reflections,
            "remove_by_index=0,3,5",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    ref = flex.reflection_table.from_file(tmp_path / "filtered.refl")
    # The test removed the reflections 0, 3 and 5
    assert len(ref) == 3
    assert list(ref["iobs"]) == [1, 2, 4]
