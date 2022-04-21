from __future__ import annotations

from dials.algorithms.profile_model.gaussian_rs.calculator import (
    _select_reflections_for_sigma_calc,
)
from dials.array_family import flex


def test_select_reflections_for_sigma_calc():
    """Test the reflection selection helper function."""
    # first select all if below threshold
    reflections = flex.reflection_table()
    reflections["id"] = flex.int(range(0, 1000))
    reflections.set_flags(flex.bool(1000, True), reflections.flags.used_in_refinement)
    reflections = _select_reflections_for_sigma_calc(reflections, 10000)
    assert reflections.size() == 1000

    # select used_in_refinement if not all used in refinement and above threshold
    reflections = flex.reflection_table()
    reflections["id"] = flex.int(range(0, 1000))
    good = flex.bool(1000, False)
    sel = flex.size_t(i for i in range(0, 1000, 2))
    good.set_selected(sel, True)
    reflections.set_flags(good, reflections.flags.used_in_refinement)
    reflections = _select_reflections_for_sigma_calc(
        reflections, min_number_of_refl=200
    )
    assert reflections.size() == 500
    assert list(reflections["id"])[0:50] == list(range(0, 100, 2))

    # top up if not enough used in refinement compared to threshold
    reflections = flex.reflection_table()
    reflections["id"] = flex.int(range(0, 1000))
    good = flex.bool(1000, False)
    sel = flex.size_t(i for i in range(0, 1000, 2))
    good.set_selected(sel, True)
    reflections.set_flags(good, reflections.flags.used_in_refinement)
    reflections = _select_reflections_for_sigma_calc(
        reflections, min_number_of_refl=700
    )
    assert reflections.size() > 700
    assert reflections.size() < 1000
