from __future__ import absolute_import, division, print_function

import procrunner
from dials.array_family import flex


def test_filter_reflections(run_in_tmpdir):
    # Make a dummy reflection table for the test setting some values and flags
    rt = flex.reflection_table.empty_standard(6)
    rt["iobs"] = flex.size_t_range(len(rt))
    rt["panel"] = flex.size_t_range(len(rt))
    rt["id"] = flex.size_t([0] * 5 + [1])
    rt["d"] = flex.double([50, 40, 3.0, 2.5, 2.0, 1.0])
    mask1 = flex.bool([True] * 3 + [False] * 3)
    mask2 = flex.bool([True, False] * 3)
    rt.set_flags(mask1, rt.flags.integrated)
    rt.set_flags(mask2, rt.flags.reference_spot)
    rt_name = "test_refs.refl"
    rt.as_pickle(rt_name)

    # Test flag expression
    cmd = [
        "dials.filter_reflections",
        rt_name,
        "flag_expression='integrated & ~reference_spot'",
    ]
    result = procrunner.run(cmd)
    assert not result.returncode and not result.stderr
    ref = flex.reflection_table.from_pickle("filtered.refl")
    # The test selects only the 2nd reflection
    assert len(ref) == 1
    assert list(ref["iobs"]) == [1]

    # Test filter by experiment id
    cmd = ["dials.filter_reflections", rt_name, "id=0"]
    result = procrunner.run(cmd)
    assert not result.returncode and not result.stderr
    ref = flex.reflection_table.from_pickle("filtered.refl")
    # The test selects only the first five reflections
    assert len(ref) == 5
    assert list(ref["iobs"]) == [0, 1, 2, 3, 4]

    # Test filter by panel
    cmd = ["dials.filter_reflections", rt_name, "panel=5"]
    result = procrunner.run(cmd)
    assert not result.returncode and not result.stderr
    ref = flex.reflection_table.from_pickle("filtered.refl")
    # The test selects only the last reflection
    assert len(ref) == 1
    assert list(ref["iobs"]) == [5]

    # Test filter by resolution
    cmd = ["dials.filter_reflections", rt_name, "d_max=3.0", "d_min=2.0"]
    result = procrunner.run(cmd)
    assert not result.returncode and not result.stderr
    ref = flex.reflection_table.from_pickle("filtered.refl")
    # The test selects only the 3rd, 4th and 5th reflections
    assert len(ref) == 3
    assert list(ref["iobs"]) == [2, 3, 4]

    # Test printing analysis
    cmd = ["dials.filter_reflections", rt_name]
    result = procrunner.run(cmd)
    assert not result.returncode and not result.stderr
