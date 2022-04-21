"""
Tests for the reflection selection algorithm.
"""

from __future__ import annotations

import itertools
from unittest.mock import Mock

from cctbx import sgtbx, uctbx
from libtbx import phil
from scitbx import sparse

from dials.algorithms.scaling.Ih_table import IhTable
from dials.algorithms.scaling.reflection_selection import (
    _loop_over_class_matrix,
    calculate_scaling_subset_ranges,
    calculate_scaling_subset_ranges_with_E2,
    select_connected_reflections_across_datasets,
)
from dials.array_family import flex


def test_select_connected_reflections_across_datasets():
    """Test the basic cross-dataset reflection selection algorithm.

    Make three reflection tables with the following reflections:
               symmetry groups
               0  1  2  3  4  5  6
            0  3  3  2  0  1  1  1
    classes 1  0  2  0  0  3  2  1
            2  2  1  1  5  0  4  0

    With target=5, expect:
    number of chosen reflections per class: [8, 7, 7]
    symmetry groups used:                   [1, 5, 0, 4]
    """

    n1 = [3, 3, 2, 0, 1, 1, 1]
    n2 = [0, 2, 0, 0, 3, 2, 1]
    n3 = [2, 1, 1, 5, 0, 4, 0]

    def make_refl_table(n_list, class_idx=0):
        """Make a reflection table with groups based on n_list."""
        r1 = flex.reflection_table()
        miller_indices = [[(0, 0, i + 1)] * n for i, n in enumerate(n_list)]
        r1["miller_index"] = flex.miller_index(
            list(itertools.chain.from_iterable(miller_indices))
        )
        r1["intensity"] = flex.double(sum(n_list), 1)
        r1["variance"] = flex.double(sum(n_list), 1)
        r1["inverse_scale_factor"] = flex.double(sum(n_list), 1)
        r1["class_index"] = flex.int(sum(n_list), class_idx)
        return r1

    reflections = [
        make_refl_table(n1, 0),
        make_refl_table(n2, 1),
        make_refl_table(n3, 2),
    ]

    space_group = sgtbx.space_group("P1")
    table = IhTable(reflections, space_group)
    experiment = Mock()
    experiment.crystal.get_space_group.return_value = space_group
    experiment.crystal.get_unit_cell.return_value = uctbx.unit_cell(
        (10, 10, 10, 90, 90, 90)
    )
    indices, datset_ids = select_connected_reflections_across_datasets(
        table, experiment, Isigma_cutoff=0.0, min_total=30, n_resolution_bins=1
    )
    assert list(indices) == [0, 1, 2, 3, 4, 5, 8, 9] + [0, 1, 2, 3, 4, 5, 6] + [
        0,
        1,
        2,
        9,
        10,
        11,
        12,
    ]
    assert list(datset_ids) == [0] * 8 + [1] * 7 + [2] * 7

    # now test again with a higher min_total
    indices, datset_ids = select_connected_reflections_across_datasets(
        table, experiment, Isigma_cutoff=0.0, min_total=10 * 3 * 4, n_resolution_bins=1
    )
    assert list(datset_ids) == [0] * 11 + [1] * 8 + [2] * 13


def test_loop_over_class_matrix():
    """Test a few different limits of the method.

    { 3, 1, 3, 2, 1, 1, 0 },
    { 2, 2, 0, 0, 3, 1, 0 },
    { 1, 4, 2, 1, 0, 0, 5 },
    """
    sorted_class_matrix = sparse.matrix(
        3,
        7,
        elements_by_columns=[
            {0: 3, 1: 2, 2: 1},
            {0: 1, 1: 2, 2: 4},
            {0: 3, 2: 2},
            {0: 2, 2: 1},
            {0: 1, 1: 3},
            {0: 1},
            {2: 5},
        ],
    )

    # first test if don't meet the minimum number
    total_in_classes, cols_not_used = _loop_over_class_matrix(
        sorted_class_matrix, 1, 8, 20
    )
    assert list(cols_not_used) == [2, 3, 4, 5, 6]
    assert list(total_in_classes) == [4.0, 4.0, 5.0]
    # now test if max per bin reached
    total_in_classes, cols_not_used = _loop_over_class_matrix(
        sorted_class_matrix, 2, 3, 8
    )
    assert list(cols_not_used) == [2, 3, 4, 5, 6]
    assert list(total_in_classes) == [4.0, 4.0, 5.0]

    # now test if request more than available - should find all
    total_in_classes, cols_not_used = _loop_over_class_matrix(
        sorted_class_matrix, 9, 50, 100
    )
    assert not cols_not_used
    assert list(total_in_classes) == [11.0, 7.0, 13.0]


def generated_param():
    """Generate a param phil scope."""
    phil_scope = phil.parse(
        """
      include scope dials.algorithms.scaling.scaling_options.phil_scope
  """,
        process_includes=True,
    )
    return phil_scope.extract()


def generated_refl_for_subset_calculation():
    """Create a reflection table suitable for splitting into blocks."""
    reflections = flex.reflection_table()
    reflections["intensity"] = flex.double([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    reflections["variance"] = flex.double(6, 1.0)
    reflections["d"] = flex.double([0.8, 2.1, 2.0, 1.4, 1.6, 2.5])
    reflections["partiality"] = flex.double(6, 1.0)
    reflections.set_flags(flex.bool(6, True), reflections.flags.integrated)
    reflections.set_flags(flex.bool(6, False), reflections.flags.bad_for_scaling)
    return reflections


def test_selection_scaling_subset_ranges_with_E2():
    """Test the scaling subset calculation with E2 range."""
    test_params = generated_param()
    rt = generated_refl_for_subset_calculation()
    rt["Esq"] = flex.double([1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    test_params.reflection_selection.E2_range = 0.8, 5.0
    test_params.reflection_selection.d_range = 1.0, 5.0  # all but first
    test_params.reflection_selection.Isigma_range = 0.9, 5.5  # all but last
    sel = calculate_scaling_subset_ranges_with_E2(rt, test_params)
    assert list(sel) == [False, True, True, True, True, False]
    rt["Esq"] = flex.double([1.0, 1.0, 1.0, 0.1, 6.0, 1.0])
    sel = calculate_scaling_subset_ranges_with_E2(rt, test_params)
    assert list(sel) == [False, True, True, False, False, False]


def test_selection_scaling_subset_ranges():
    """Test the scaling subset calculation with E2 range."""
    test_params = generated_param()
    rt = generated_refl_for_subset_calculation()
    test_params.reflection_selection.E2_range = 0.8, 5.0
    test_params.reflection_selection.d_range = 1.0, 5.0  # all but first
    test_params.reflection_selection.Isigma_range = 0.9, 5.5  # all but last
    sel = calculate_scaling_subset_ranges(rt, test_params)
    assert list(sel) == [False, True, True, True, True, False]
