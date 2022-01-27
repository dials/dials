from __future__ import annotations

from cctbx import sgtbx, uctbx
from scitbx.array_family import flex

from dials.algorithms.integration import filtering


def test_powder_ring_filter():
    unit_cell = uctbx.unit_cell((4.498, 4.498, 7.338, 90, 90, 120))
    space_group = sgtbx.space_group_info(number=194).group()
    d_spacings = flex.double([1.9220704466392748])
    d_min = flex.min(d_spacings)
    ice_filter = filtering.PowderRingFilter(unit_cell, space_group, d_min, width=0.004)
    assert min(uctbx.d_star_sq_as_d(ice_filter.d_star_sq)) < d_min
    assert all(ice_filter(d_spacings))
