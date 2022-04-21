"""
tests for functions in dials.util.exclude_images.py
"""

from __future__ import annotations

import copy
from unittest.mock import Mock

import pytest

from dxtbx.model import Experiment, ExperimentList, Scan

from dials.array_family import flex
from dials.util.exclude_images import (
    _parse_exclude_images_commands,
    exclude_image_ranges_for_scaling,
    exclude_image_ranges_from_scans,
    get_selection_for_valid_image_ranges,
    get_valid_image_ranges,
    set_initial_valid_image_ranges,
)


def make_scan_experiment(image_range=(1, 100), expid="0"):
    """Make an experiment with a scan"""
    return Experiment(scan=Scan(image_range, (0.0, 1.0)), identifier=expid)


def make_scanless_experiment(expid="1"):
    """Make an experiment without a scan i.e. single-image dataset"""
    return Experiment(identifier=expid)


def test_parse_exclude_images_commands():
    """Test for namesake function"""
    formats = (
        [["1:101:200"], ["0:201:300"]],  # if given as separate exclude_images=
        [["1:101:200,0:201:300"]],  # if given as exclude_images="1:101:200,0:201:300"
        [
            ["1:101:200", "0:201:300"]
        ],  # if given as exclude_images="1:101:200 0:201:300"
    )
    for command in formats:
        r1 = flex.reflection_table()
        r1.experiment_identifiers()[1] = "1"
        r0 = flex.reflection_table()
        r0.experiment_identifiers()[0] = "0"
        tables = [r0, r1]
        ranges = _parse_exclude_images_commands(command, [], tables)
        assert ranges == [("1", (101, 200)), ("0", (201, 300))]

    experiments = ["1", "2"]
    short_command = [["101:200"]]
    with pytest.raises(ValueError):
        _ = _parse_exclude_images_commands(short_command, experiments, tables)
    mock_exp = Mock()
    mock_exp.identifier = "1"
    ranges = _parse_exclude_images_commands(short_command, [mock_exp], tables)
    assert ranges == [("1", (101, 200))]
    with pytest.raises(ValueError):
        _ = _parse_exclude_images_commands([["1:101-200"]], [mock_exp], tables)
    with pytest.raises(ValueError):
        _ = _parse_exclude_images_commands([["1:101:a"]], [], tables)


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
    r1 = flex.reflection_table()
    r1.experiment_identifiers()[1] = "1"
    r0 = flex.reflection_table()
    r0.experiment_identifiers()[0] = "0"
    tables = [r0, r1]
    explist = exclude_image_ranges_from_scans(tables, explist, exclude_images)
    assert list(explist[0].scan.get_valid_image_ranges("0")) == [(1, 80)]
    assert list(explist[1].scan.get_valid_image_ranges("1")) == [(1, 60), (81, 100)]
    # Try excluding a range that already has been excluded
    explist = exclude_image_ranges_from_scans(tables, explist, [["1:70:80"]])
    assert list(explist[0].scan.get_valid_image_ranges("0")) == [(1, 80)]
    assert list(explist[1].scan.get_valid_image_ranges("1")) == [(1, 60), (81, 100)]
    scanlessexplist = ExperimentList([make_scanless_experiment()])
    with pytest.raises(ValueError):
        _ = exclude_image_ranges_from_scans(tables, scanlessexplist, [["0:1:100"]])
    # Now try excluding everything, should set an empty array
    explist = exclude_image_ranges_from_scans(tables, explist, [["1:1:100"]])
    assert list(explist[0].scan.get_valid_image_ranges("0")) == [(1, 80)]
    assert list(explist[1].scan.get_valid_image_ranges("1")) == []

    ## test what happens if a single image is left within the scan
    explist = ExperimentList(
        [make_scan_experiment(expid="0"), make_scan_experiment(expid="1")]
    )
    exclude_images = [["0:81:100"], ["1:76:79"], ["1:81:99"]]
    r1 = flex.reflection_table()
    r1.experiment_identifiers()[1] = "1"
    r0 = flex.reflection_table()
    r0.experiment_identifiers()[0] = "0"
    tables = [r0, r1]
    explist = exclude_image_ranges_from_scans(tables, explist, exclude_images)
    assert list(explist[0].scan.get_valid_image_ranges("0")) == [(1, 80)]
    assert list(explist[1].scan.get_valid_image_ranges("1")) == [
        (1, 75),
        (80, 80),
        (100, 100),
    ]


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
    refl1.experiment_identifiers()[0] = "0"
    refl2.experiment_identifiers()[1] = "1"
    explist = ExperimentList(
        [
            make_scan_experiment(image_range=(2, 20), expid="0"),
            make_scan_experiment(image_range=(2, 20), expid="1"),
        ]
    )
    refls, explist = exclude_image_ranges_for_scaling(
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
