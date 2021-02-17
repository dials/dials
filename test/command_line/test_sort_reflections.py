import procrunner


def test_sort_intensities(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dev.dials.sort_reflections",
            dials_data("centroid_test_data").join("integrated.pickle").strpath,
            "key=intensity.sum.value",
            "output=sorted1.refl",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("sorted1.refl").check(file=1)

    from dials.array_family import flex

    data = flex.reflection_table.from_file(tmpdir.join("sorted1.refl").strpath)
    assert_sorted(data["intensity.sum.value"])


def test_reverse_sort_intensities(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dev.dials.sort_reflections",
            dials_data("centroid_test_data").join("integrated.pickle").strpath,
            "output=sorted2.refl",
            "key=intensity.sum.value",
            "reverse=True",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("sorted2.refl").check(file=1)

    from dials.array_family import flex

    data = flex.reflection_table.from_file(tmpdir.join("sorted2.refl").strpath)
    assert_sorted(data["intensity.sum.value"], reverse=True)


def test_default_sort_on_miller_index(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dev.dials.sort_reflections",
            dials_data("centroid_test_data").join("integrated.pickle").strpath,
            "output=sorted3.refl",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("sorted3.refl").check(file=1)

    from dials.array_family import flex

    data = flex.reflection_table.from_file(tmpdir.join("sorted3.refl").strpath)
    mi1 = data["miller_index"]
    orig = flex.reflection_table.from_file(
        dials_data("centroid_test_data").join("integrated.pickle").strpath
    )
    mi2 = flex.miller_index(sorted(orig["miller_index"]))
    assert mi1.all_eq(mi2)


def test_default_sort_on_miller_index_verbose(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dev.dials.sort_reflections",
            dials_data("centroid_test_data").join("integrated.pickle").strpath,
            "output=sorted4.refl",
            "-v",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("sorted4.refl").check(file=1)
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
