"""Tests for the BinnedStatistics class"""

from __future__ import annotations

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
    assert binned_statistics.bin_is_sorted().all_eq(
        flex.bool((False, False, False, False))
    )


def test_median_and_iqr():

    # Even number of values, one bin
    vals1 = flex.double((1, 2, 3, 4))
    bins = flex.size_t([0] * len(vals1))
    binned_statistics = BinnedStatistics(vals1, bins, 1)
    assert binned_statistics.get_medians().all_eq(flex.double((2.5,)))
    assert binned_statistics.get_iqrs().all_eq(flex.double((1.5,)))
    # NB this differs from the calculation performed via
    # _,q1,_,q3,_=five_number_summary(vals1);q3-q1
    # which gives 2.0. However, it matches the result of R's IQR function

    # Odd number of values, one bin, and robustness test
    vals2 = flex.double((1, 2, 3, 100, 1000))
    bins = flex.size_t([0] * len(vals2))
    binned_statistics = BinnedStatistics(vals2, bins, 1)
    assert binned_statistics.get_medians().all_eq(flex.double((3.0,)))
    assert binned_statistics.get_iqrs().all_eq(flex.double((98,)))
    # NB this is the same as the calculation performed via
    # _,q1,_,q3,_=five_number_summary(vals2);q3-q1, and matches R's IQR function

    # Now combine the data and randomise order
    vals3 = vals1.concatenate(vals2)
    bins = flex.size_t([0] * len(vals1) + [1] * len(vals2))
    perm = flex.random_permutation(len(bins))
    vals3 = vals3.select(perm)
    bins = bins.select(perm)
    binned_statistics = BinnedStatistics(vals3, bins, 2)
    assert binned_statistics.get_medians().all_eq(
        flex.double(
            (
                2.5,
                3.0,
            )
        )
    )
    assert binned_statistics.get_iqrs().all_eq(
        flex.double(
            (
                1.5,
                98,
            )
        )
    )
