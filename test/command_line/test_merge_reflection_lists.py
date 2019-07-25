from __future__ import absolute_import, division, print_function

import procrunner


def test(run_in_tmpdir):
    from dials.array_family import flex

    table = flex.reflection_table()
    table["hkl"] = flex.miller_index(360)
    table["id"] = flex.int(360)
    table["intensity.sum.value"] = flex.double(360)
    table.as_pickle("temp1.refl")
    table.as_pickle("temp2.refl")

    result = procrunner.run(
        [
            "dev.dials.merge_reflection_lists",
            "temp1.refl",
            "temp2.refl",
            "method=update",
        ]
    )
    assert not result.returncode and not result.stderr

    table = flex.reflection_table.from_pickle("merged.refl")
    assert len(table) == 360

    result = procrunner.run(
        [
            "dev.dials.merge_reflection_lists",
            "temp1.refl",
            "temp2.refl",
            "method=extend",
        ]
    )
    assert not result.returncode and not result.stderr

    table = flex.reflection_table.from_pickle("merged.refl")
    assert len(table) == 720
