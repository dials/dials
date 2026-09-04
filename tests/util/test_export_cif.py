"""Tests for the parts of the small molecule CIF exporter that are awkward to
reach through dials.export itself."""

from __future__ import annotations

import pytest

from cctbx import sgtbx, uctbx

from dials.util.export_cif import _constrained_angles, _theta_from_d


@pytest.mark.parametrize(
    "space_group,cell,expected",
    [
        ("P 1", (5, 8, 12, 88, 92, 105), [False, False, False]),
        ("P 1 21 1", (5, 8, 12, 90, 92, 90), [True, False, True]),
        ("P 2 2 2", (5, 8, 12, 90, 90, 90), [True, True, True]),
        ("P 4", (5, 5, 12, 90, 90, 90), [True, True, True]),
        ("R 3 :H", (5, 5, 12, 90, 90, 120), [True, True, True]),
    ],
)
def test_constrained_angles(space_group, cell, expected):
    sg = sgtbx.space_group_info(space_group).group()
    assert _constrained_angles(sg, uctbx.unit_cell(cell)) == expected


def test_theta_from_d():
    # 2 d sin(theta) = lambda
    assert _theta_from_d(1.0, 1.0) == pytest.approx(30.0)
    assert _theta_from_d(10.0, 1.0) == pytest.approx(2.8660, abs=1e-4)
    # a d-spacing below lambda / 2 is not accessible, so clamp at 90 degrees
    assert _theta_from_d(0.1, 1.0) == pytest.approx(90.0)
