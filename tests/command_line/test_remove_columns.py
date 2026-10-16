from __future__ import annotations

import shutil
import subprocess

import pytest

from dxtbx.util import ersatz_uuid4

from dials.command_line.remove_columns import expand_column_patterns, format_file_size

COLUMNS = [
    "hkl",
    "id",
    "panel",
    "xyzobs.px.value",
    "xyzobs.mm.value",
    "xyzobs.mm.variance",
    "intensity.sum.value",
]


@pytest.mark.parametrize(
    "commands,expected,unmatched",
    [
        (
            [["xyzobs.mm.value,xyzobs.mm.variance"]],
            {"xyzobs.mm.value", "xyzobs.mm.variance"},
            [],
        ),
        (
            [["xyzobs.mm.value", "xyzobs.mm.variance"]],
            {"xyzobs.mm.value", "xyzobs.mm.variance"},
            [],
        ),
        (
            [["xyzobs.mm.value"], ["xyzobs.mm.variance"]],
            {"xyzobs.mm.value", "xyzobs.mm.variance"},
            [],
        ),
        ([["xyzobs.mm.*"]], {"xyzobs.mm.value", "xyzobs.mm.variance"}, []),
        ([["hkl,nope"]], {"hkl"}, ["nope"]),
        ([None], set(), []),
        ([], set(), []),
    ],
)
def test_expand_column_patterns(commands, expected, unmatched):
    assert expand_column_patterns(commands, COLUMNS) == (expected, unmatched)


@pytest.mark.parametrize(
    "nbytes,expected",
    [(0, "0 B"), (512, "512 B"), (2048, "2.0 kB"), (1234567, "1.2 MB")],
)
def test_format_file_size(nbytes, expected):
    assert format_file_size(nbytes) == expected


@pytest.fixture
def reflection_file(tmp_path):
    """A small reflection table with a representative set of columns."""

    from dials.array_family import flex

    table = flex.reflection_table()
    table["hkl"] = flex.miller_index(100)
    table["id"] = flex.int(100)
    table["panel"] = flex.size_t(100)
    table["xyzobs.px.value"] = flex.vec3_double(100)
    table["xyzobs.mm.value"] = flex.vec3_double(100)
    table["xyzobs.mm.variance"] = flex.vec3_double(100)
    table["intensity.sum.value"] = flex.double(100)
    for identifier in set(table["id"]):
        table.experiment_identifiers()[identifier] = ersatz_uuid4()
    table.as_file(tmp_path / "input.refl")

    return tmp_path / "input.refl"


def test_remove_named_columns(tmp_path, reflection_file):
    from dials.array_family import flex

    result = subprocess.run(
        [
            shutil.which("dials.remove_columns"),
            reflection_file,
            "remove=xyzobs.mm.value,xyzobs.mm.variance",
            "output=small.refl",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "small.refl").is_file()

    table = flex.reflection_table.from_file(tmp_path / "small.refl")
    assert len(table) == 100
    assert "xyzobs.mm.value" not in table
    assert "xyzobs.mm.variance" not in table
    assert set(table.keys()) == {
        "hkl",
        "id",
        "panel",
        "xyzobs.px.value",
        "intensity.sum.value",
    }
    # experiment identifiers should survive
    assert list(table.experiment_identifiers().keys()) == [0]


def test_default_output_filename(tmp_path, reflection_file):
    from dials.array_family import flex

    result = subprocess.run(
        [
            shutil.which("dials.remove_columns"),
            reflection_file,
            "remove=intensity.sum.value",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    table = flex.reflection_table.from_file(tmp_path / "stripped.refl")
    assert "intensity.sum.value" not in table
    assert "hkl" in table


def test_remove_columns_wildcard(tmp_path, reflection_file):
    from dials.array_family import flex

    result = subprocess.run(
        [
            shutil.which("dials.remove_columns"),
            reflection_file,
            "remove=xyzobs.mm.*",
            "output=small.refl",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    table = flex.reflection_table.from_file(tmp_path / "small.refl")
    assert "xyzobs.mm.value" not in table
    assert "xyzobs.mm.variance" not in table
    # the px column matches neither pattern and must be retained
    assert "xyzobs.px.value" in table


def test_remove_columns_repeated_option(tmp_path, reflection_file):
    from dials.array_family import flex

    result = subprocess.run(
        [
            shutil.which("dials.remove_columns"),
            reflection_file,
            "remove=xyzobs.mm.value",
            "remove=xyzobs.mm.variance",
            "output=small.refl",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    table = flex.reflection_table.from_file(tmp_path / "small.refl")
    assert "xyzobs.mm.value" not in table
    assert "xyzobs.mm.variance" not in table


def test_unknown_column_warns_but_succeeds(tmp_path, reflection_file):
    from dials.array_family import flex

    result = subprocess.run(
        [
            shutil.which("dials.remove_columns"),
            reflection_file,
            "remove=not_a_column,intensity.sum.value",
            "output=small.refl",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert b"not_a_column" in result.stdout

    table = flex.reflection_table.from_file(tmp_path / "small.refl")
    assert "intensity.sum.value" not in table


def test_no_columns_matched_is_an_error(tmp_path, reflection_file):
    result = subprocess.run(
        [
            shutil.which("dials.remove_columns"),
            reflection_file,
            "remove=not_a_column",
            "output=small.refl",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert result.returncode
    assert not (tmp_path / "small.refl").is_file()


def test_no_remove_option_is_an_error(tmp_path, reflection_file):
    result = subprocess.run(
        [shutil.which("dials.remove_columns"), reflection_file],
        cwd=tmp_path,
        capture_output=True,
    )
    assert result.returncode
    assert not (tmp_path / "stripped.refl").is_file()


def test_removing_every_column_is_an_error(tmp_path, reflection_file):
    result = subprocess.run(
        [
            shutil.which("dials.remove_columns"),
            reflection_file,
            "remove=*",
            "output=small.refl",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert result.returncode
    assert not (tmp_path / "small.refl").is_file()
