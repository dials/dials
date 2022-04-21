"""
Test for the functions in dials.util.batch_handling.

These functions are designed to handle reflection/experiment lists that
include both sequence experiments and single image datasets, which do not
have a scan object.
"""

from __future__ import annotations

from dxtbx.model import Experiment, Scan

from dials.array_family import flex
from dials.util.batch_handling import (
    _next_epoch,
    assign_batches_to_reflections,
    calculate_batch_offsets,
    get_batch_ranges,
    get_image_ranges,
    set_batch_offsets,
)


def reflections_1():
    """Test reflection table with batch"""
    r = flex.reflection_table()
    r["batch"] = flex.int([1, 11, 21, 31, 41, 51, 61, 71, 81, 91])
    r.set_flags(flex.bool(10, False), r.flags.user_excluded_in_scaling)
    return r


def reflections_2():
    """Test reflection table with batch"""
    r = flex.reflection_table()
    r["batch"] = flex.int([201, 211, 221, 231, 241, 251, 261, 271, 281, 291])
    r.set_flags(flex.bool(10, False), r.flags.user_excluded_in_scaling)
    return r


def reflections_3():
    """Test reflection table with xyzobs.px.value"""
    r = flex.reflection_table()
    r["xyzobs.px.value"] = flex.vec3_double([(0, 0, 0.5), (0, 0, 1.5)])
    return r


def test_assign_batches_to_reflections():
    """Test for namesake function"""
    reflections = [reflections_3(), reflections_3()]
    reflections = assign_batches_to_reflections(reflections, batch_offsets=[0, 100])
    assert list(reflections[0]["batch"]) == [1, 2]
    assert list(reflections[1]["batch"]) == [101, 102]


def test_calculate_batch_offsets():
    """Test offset calculation. Offset is next number ending in 01 bigger than
    previous batch numbers which is not consecutive"""
    scan = Scan(image_range=[1, 200], oscillation=[0.0, 1.0])
    exp1 = Experiment(scan=scan)
    exp2 = Experiment()
    offsets = calculate_batch_offsets([exp1, exp2])
    assert offsets == [0, 301]


def test_set_batch_offsets():
    """Test for namesake function"""
    scan = Scan(image_range=[1, 200], oscillation=[0.0, 1.0])
    exp1 = Experiment(scan=scan)
    exp2 = Experiment()
    set_batch_offsets([exp2, exp1], [0, 101])
    assert exp1.scan.get_batch_offset() == 101


def test_get_batch_ranges():
    """Test for namesake function"""
    scan = Scan(image_range=[1, 200], oscillation=[0.0, 1.0])
    exp1 = Experiment(scan=scan)
    exp2 = Experiment(scan=scan)
    batch_offsets = [0, 300]
    experiments = [exp1, exp2]
    batch_ranges = get_batch_ranges(experiments, batch_offsets)
    assert batch_ranges == [(1, 200), (301, 500)]


def test_get_image_ranges():
    """Test for namesake function"""
    scan = Scan(image_range=[1, 200], oscillation=[0.0, 1.0])
    exp1 = Experiment(scan=scan)
    exp2 = Experiment()
    experiments = [exp1, exp2]
    image_ranges = get_image_ranges(experiments)
    assert image_ranges == [(1, 200), (0, 0)]


def test_next_epoch():
    """Test for namesake function"""
    assert _next_epoch(100) == 201
    assert _next_epoch(99) == 101
    assert _next_epoch(105) == 201
