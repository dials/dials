from __future__ import absolute_import, division, print_function

import procrunner


def test(run_in_tmpdir):
    from dials.array_family import flex

    table = flex.reflection_table()
    table["hkl"] = flex.miller_index(360)
    table["id"] = flex.int(360)
    table["intensity.sum.value"] = flex.double(360)
    table.as_file("temp1.pickle")
    table.as_file("temp2.pickle")

    result = procrunner.run(
        [
            "dev.dials.merge_reflection_lists",
            "temp1.pickle",
            "temp2.pickle",
            "method=update",
        ]
    )
    assert result["exitcode"] == 0
    assert result["stderr"] == ""

    table = flex.reflection_table.from_file("merged.pickle")
    assert len(table) == 360

    result = procrunner.run(
        [
            "dev.dials.merge_reflection_lists",
            "temp1.pickle",
            "temp2.pickle",
            "method=extend",
        ]
    )
    assert result["exitcode"] == 0
    assert result["stderr"] == ""

    table = flex.reflection_table.from_file("merged.pickle")
    assert len(table) == 720
