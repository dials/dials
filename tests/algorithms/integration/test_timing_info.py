from __future__ import annotations

from dials.algorithms.integration import TimingInfo


def test_timing_info():
    """Create TimingInfo objects, test addition and str methods."""
    time1 = TimingInfo()
    time2 = TimingInfo()

    time1.read = 5
    time1.process = 4
    time1.user = 10
    time2.read = 15
    time2.process = 6
    time2.user = 1

    total = time1 + time2
    assert total.read == 20
    assert total.process == 10
    assert total.user == 11

    s = str(total)
    assert (
        s
        == """+--------------+---------------+
| Read time    | 20.00 seconds |
| Process time | 10.00 seconds |
| User time    | 11.00 seconds |
+--------------+---------------+"""
    )
