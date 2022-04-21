from __future__ import annotations

import procrunner

from dials.array_family import flex


def test_sort_intensities(dials_data, tmp_path):
    result = procrunner.run(
        [
            "dials.sort_reflections",
            dials_data("centroid_test_data", pathlib=True) / "integrated.pickle",
            "key=intensity.sum.value",
            "output=sorted1.refl",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("sorted1.refl").is_file()

    data = flex.reflection_table.from_file(tmp_path / "sorted1.refl")
    assert_sorted(data["intensity.sum.value"])


def test_reverse_sort_intensities(dials_data, tmp_path):
    result = procrunner.run(
        [
            "dials.sort_reflections",
            dials_data("centroid_test_data", pathlib=True) / "integrated.pickle",
            "output=sorted2.refl",
            "key=intensity.sum.value",
            "reverse=True",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("sorted2.refl").is_file()

    data = flex.reflection_table.from_file(tmp_path / "sorted2.refl")
    assert_sorted(data["intensity.sum.value"], reverse=True)


def test_default_sort_on_miller_index(dials_data, tmp_path):
    result = procrunner.run(
        [
            "dials.sort_reflections",
            dials_data("centroid_test_data", pathlib=True) / "integrated.pickle",
            "output=sorted3.refl",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("sorted3.refl").is_file()

    data = flex.reflection_table.from_file(tmp_path / "sorted3.refl")
    mi1 = data["miller_index"]
    orig = flex.reflection_table.from_file(
        dials_data("centroid_test_data", pathlib=True) / "integrated.pickle"
    )
    mi2 = flex.miller_index(sorted(orig["miller_index"]))
    assert mi1.all_eq(mi2)


def test_default_sort_on_miller_index_verbose(dials_data, tmp_path):
    result = procrunner.run(
        [
            "dials.sort_reflections",
            dials_data("centroid_test_data", pathlib=True) / "integrated.pickle",
            "output=sorted4.refl",
            "-v",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("sorted4.refl").is_file()
    assert b"Head of sorted list miller_index:" in result.stdout


def assert_sorted(data, reverse=False):
    assert len(data) > 0
    elem = data[0]
    for x in data:
        if reverse is True:
            assert x <= elem
        else:
            assert x >= elem
        elem = x
