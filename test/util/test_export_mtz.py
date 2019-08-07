"""
Unit testing for the export_mtz.py routines
"""

from __future__ import absolute_import, division, print_function

import itertools

from dials.util.batch_handling import calculate_batch_offsets


def in_ranges(value, ranges):
    """Check if a value is in a list of ranges"""
    return all(low <= value <= high for low, high in ranges)


def range_to_set(ranges):
    """Convert a list of ranges to a set of numbers"""
    return set().union(*[set(range(l, h + 1)) for l, h in ranges])


def has_consecutive_ranges(ranges):
    """Check that a set of ranges is non-consecutive"""
    total_entries = sum(h - l + 1 + 1 for l, h in ranges)
    union = set().union(*[set(range(l, h + 2)) for l, h in ranges])
    return len(union) < total_entries


def offset_ranges(offsets, ranges):
    return [(l + off, h + off) for off, (l, h) in zip(offsets, ranges)]


class TestBatchRangeCalculations(object):
    "Test the calculation of non-overlapping batch ranges"

    class MockScan(object):
        def __init__(self, image_range):
            self._batch_offset = 0
            self._image_range = image_range

        def get_batch_offset(self):
            return self._batch_offset

        def set_batch_offset(self, batch_offset):
            self._batch_offset = batch_offset

        def get_image_range(self):
            return self._image_range

    class MockExperiment(object):
        def __init__(self, image_range, scan=True):
            assert len(image_range) == 2
            self.scaling_model = None
            if scan:
                self.scan = TestBatchRangeCalculations.MockScan(image_range)
            else:
                self.scan = []

    def _run_ranges(self, ranges):
        """Convenience method to run the routine with a minimal experiment, and return the result as ranges of batch number"""
        input_data = [self.MockExperiment(x) for x in ranges]
        return offset_ranges(calculate_batch_offsets(input_data), ranges)

    def _run_ranges_to_set(self, ranges):
        """Runs a list of ranges and returns a set of individual batch numbers"""
        return range_to_set(self._run_ranges(ranges))

    def test_calculate_batch_ranges(self):
        assert self._run_ranges([(1, 1)]) == [(1, 1)]

        # Zero is shifted
        assert all(
            [x > 0 for x in self._run_ranges_to_set([(0, 0)])]
        ), "Should be no zeroth/negative batch"

        assert not set(self._run_ranges([(1, 1), (1, 1)])) == {
            (1, 1)
        }, "Overlapping simple ranges"

        data_tests = [
            [(1, 1), (1, 1)],
            # while we decide behaviour, remove input
            #            [(1, 1), (8, 8), (9, 9)],
            [(23, 24), (70, 100), (1, 1), (1, 4), (1, 1)],
            [(0, 98)],
        ]
        for data in data_tests:
            print("Running ", data)
            print("  ", self._run_ranges(data))
            assert all(
                [
                    float(x).is_integer()
                    for x in itertools.chain(*self._run_ranges(data))
                ]
            ), "Fractional epochs"
            assert all(
                isinstance(x, int) for x in itertools.chain(*self._run_ranges(data))
            ), "Not all true integers"
            assert all(
                [x > 0 for x in self._run_ranges_to_set([(0, 0)])]
            ), "Should be no zeroth/negative batch"
            assert not has_consecutive_ranges(self._run_ranges(data))

        exp1 = TestBatchRangeCalculations.MockExperiment((1, 1), scan=False)
        exp2 = TestBatchRangeCalculations.MockExperiment((1, 1), scan=False)
        offsets = calculate_batch_offsets([exp1, exp2])
        assert all(float(x).is_integer() for x in offsets)
        assert all(isinstance(x, int) for x in offsets)
        assert all(x > 0 for x in offsets)
