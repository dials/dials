"""Tests for the BinnedStatistics class"""

from scitbx.array_family import flex

from dials.algorithms.statistics import BinnedStatistics


def test_partitions():

    vals = flex.double((1, 2, 3, 1, 2, 3, 1, 2, 3))
    bins = flex.size_t((0, 1, 3, 0, 1, 3, 0, 1, 3))
    binned_statistics = BinnedStatistics(vals, bins, 4)

    assert binned_statistics.get_values_in_bin(0).all_eq(flex.double((1, 1, 1)))
    assert binned_statistics.get_values_in_bin(1).all_eq(flex.double((2, 2, 2)))
    assert binned_statistics.get_values_in_bin(2).size() == 0
    assert binned_statistics.get_values_in_bin(3).all_eq(flex.double((3, 3, 3)))
    assert binned_statistics.bin_is_empty().all_eq(
        flex.bool((False, False, True, False))
    )
