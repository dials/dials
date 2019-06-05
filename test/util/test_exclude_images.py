"""
tests for functions in dials.util.exclude_images.py
"""
from __future__ import absolute_import, division, print_function
import copy
import pytest
from dxtbx.model import Experiment, Scan, ExperimentList
from dials.array_family import flex
from dials.util.exclude_images import (
    _parse_exclude_images_commands,
    set_initial_valid_image_ranges,
    get_valid_image_ranges,
    exclude_image_ranges_from_scans,
    get_selection_for_valid_image_ranges,
    exclude_image_ranges_for_scaling,
)


def make_scan_experiment(image_range=(1, 100), expid="0"):
    """Make an experiment with a scan"""
    return Experiment(scan=Scan(image_range, (0.0, 1.0)), identifier=expid)


def make_scanless_experiment(expid="1"):
    """Make an experiment without a scan i.e. single-image dataset"""
    return Experiment(identifier=expid)


def test_parse_exclude_images_commands():
    """Test for namesake function"""
    commands = [["1:101:200"], ["0:201:300"]]
    ranges = _parse_exclude_images_commands(commands)
    assert ranges == [("1", (101, 200)), ("0", (201, 300))]
    with pytest.raises(ValueError):
        _ = _parse_exclude_images_commands([["1:101-200"]])
    with pytest.raises(ValueError):
        _ = _parse_exclude_images_commands([["1:101"]])
    with pytest.raises(ValueError):
        _ = _parse_exclude_images_commands([["1:101:a"]])


def test_set_get_initial_valid_image_ranges():
    """Test for get/set valid_image_ranges functions"""
    explist = ExperimentList([make_scan_experiment(), make_scanless_experiment()])
    explist = set_initial_valid_image_ranges(explist)
    assert list(explist[0].scan.get_valid_image_ranges("0")) == [(1, 100)]
    ranges = get_valid_image_ranges(explist)
    assert len(ranges) == 2
    assert list(ranges[0]) == [(1, 100)]
    assert ranges[1] is None


def test_exclude_image_ranges_from_scans():
    """Test for namesake function"""
    explist = ExperimentList(
        [make_scan_experiment(expid="0"), make_scan_experiment(expid="1")]
    )
    exclude_images = [["0:81:100"], ["1:61:80"]]
    explist = exclude_image_ranges_from_scans(explist, exclude_images)
    assert list(explist[0].scan.get_valid_image_ranges("0")) == [(1, 80)]
    assert list(explist[1].scan.get_valid_image_ranges("1")) == [(1, 60), (81, 100)]
    # Try excluding a range that already has been excluded
    explist = exclude_image_ranges_from_scans(explist, [["1:70:80"]])
    assert list(explist[0].scan.get_valid_image_ranges("0")) == [(1, 80)]
    assert list(explist[1].scan.get_valid_image_ranges("1")) == [(1, 60), (81, 100)]
    scanlessexplist = ExperimentList([make_scanless_experiment()])
    with pytest.raises(ValueError):
        _ = exclude_image_ranges_from_scans(scanlessexplist, [["0:1:100"]])
    # Now try excluding everything, should set an empty array
    explist = exclude_image_ranges_from_scans(explist, [["1:1:100"]])
    assert list(explist[0].scan.get_valid_image_ranges("0")) == [(1, 80)]
    assert list(explist[1].scan.get_valid_image_ranges("1")) == []


def test_get_selection_for_valid_image_ranges():
    """Test for namesake function"""
    exp = make_scan_experiment()
    exp.scan.set_valid_image_ranges("0", [(2, 10)])
    refl = flex.reflection_table()
    refl["xyzobs.px.value"] = flex.vec3_double(
        [(0, 0, 0.5), (0, 0, 1.5), (0, 0, 5.5), (0, 0, 9.5), (0, 0, 10.5)]
    )
    sel = get_selection_for_valid_image_ranges(refl, exp)
    assert list(sel) == [False, True, True, True, False]
    exp = make_scanless_experiment()
    assert list(get_selection_for_valid_image_ranges(refl, exp)) == [True] * 5


def test_exclude_image_ranges_for_scaling():
    """Test for namesake function."""
    refl1 = flex.reflection_table()
    refl1["xyzobs.px.value"] = flex.vec3_double(
        [(0, 0, 0.5), (0, 0, 1.5), (0, 0, 5.5), (0, 0, 9.5), (0, 0, 10.5)]
    )
    refl1.set_flags(flex.bool(5, False), refl1.flags.user_excluded_in_scaling)
    refl2 = copy.deepcopy(refl1)
    explist = ExperimentList(
        [
            make_scan_experiment(image_range=(2, 20), expid="0"),
            make_scan_experiment(image_range=(2, 20), expid="1"),
        ]
    )
    refls, exps = exclude_image_ranges_for_scaling(
        [refl1, refl2], explist, [["1:11:20"]]
    )
    assert list(explist[0].scan.get_valid_image_ranges("0")) == [(2, 20)]
    assert list(explist[1].scan.get_valid_image_ranges("1")) == [(2, 10)]
    assert list(refls[0].get_flags(refls[0].flags.user_excluded_in_scaling)) == [
        True,
        False,
        False,
        False,
        False,
    ]
    assert list(refls[1].get_flags(refls[0].flags.user_excluded_in_scaling)) == [
        True,
        False,
        False,
        False,
        True,
    ]
