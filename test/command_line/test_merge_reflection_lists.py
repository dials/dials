from __future__ import absolute_import, division, print_function

import procrunner


def test(run_in_tmpdir):
    from dials.array_family import flex

    table = flex.reflection_table()
    table["hkl"] = flex.miller_index(360)
    table["id"] = flex.int(360)
    table["intensity.sum.value"] = flex.double(360)
    table.as_msgpack_file("temp1.mpack")
    table.as_msgpack_file("temp2.mpack")

    result = procrunner.run(
        [
            "dev.dials.merge_reflection_lists",
            "temp1.mpack",
            "temp2.mpack",
            "method=update",
        ]
    )
    assert result["exitcode"] == 0
    assert result["stderr"] == ""

    table = flex.reflection_table.from_msgpack_file("merged.mpack")
    assert len(table) == 360

    result = procrunner.run(
        [
            "dev.dials.merge_reflection_lists",
            "temp1.mpack",
            "temp2.mpack",
            "method=extend",
        ]
    )
    assert result["exitcode"] == 0
    assert result["stderr"] == ""

    table = flex.reflection_table.from_msgpack_file("merged.mpack")
    assert len(table) == 720
