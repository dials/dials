from __future__ import annotations

import pickle
import random

import pytest

from dials.algorithms.background.modeller import (
    BackgroundStatistics,
    MultiPanelBackgroundStatistics,
)
from dials.array_family import flex
from dials.model.data import ImageVolume, MultiPanelImageVolume


def make_image_volume(frame0, frame1, height, width):
    volume = ImageVolume(frame0, frame1, height, width)
    for frame in range(frame0, frame1):
        data = flex.double([random.uniform(0, 100) for _ in range(height * width)])
        data.reshape(flex.grid(height, width))
        mask = flex.bool([random.random() < 0.8 for _ in range(height * width)])
        mask.reshape(flex.grid(height, width))
        volume.set_image(frame, data, mask)
    return volume


def assert_statistics_equal(a, b):
    for array in ("sum", "sum_sq", "num", "min", "max"):
        x, y = getattr(a, array)(), getattr(b, array)()
        assert y.all() == x.all()
        assert list(y) == list(x)
    for array in ("mean", "variance", "dispersion", "mask"):
        x, y = getattr(a, array)(1), getattr(b, array)(1)
        assert list(y) == list(x)


@pytest.mark.parametrize("protocol", range(pickle.HIGHEST_PROTOCOL + 1))
def test_background_statistics_pickle(protocol):
    stats = BackgroundStatistics(make_image_volume(0, 5, 6, 7))
    stats2 = pickle.loads(pickle.dumps(stats, protocol=protocol))
    assert_statistics_equal(stats, stats2)


def test_unpickled_background_statistics_can_be_accumulated():
    volume1 = make_image_volume(0, 3, 4, 5)
    volume2 = make_image_volume(3, 6, 4, 5)

    expected = BackgroundStatistics(volume1)
    expected += BackgroundStatistics(volume2)

    # As done when combining results from worker processes
    result = pickle.loads(pickle.dumps(BackgroundStatistics(volume1)))
    result += pickle.loads(pickle.dumps(BackgroundStatistics(volume2)))
    assert_statistics_equal(expected, result)


def test_background_statistics_from_arrays():
    stats = BackgroundStatistics(make_image_volume(0, 3, 4, 5))
    stats2 = BackgroundStatistics(
        sum=stats.sum(),
        sum_sq=stats.sum_sq(),
        num=stats.num(),
        min=stats.min(),
        max=stats.max(),
    )
    assert_statistics_equal(stats, stats2)

    other = BackgroundStatistics(make_image_volume(0, 3, 5, 4))
    with pytest.raises(RuntimeError):
        BackgroundStatistics(
            stats.sum(), stats.sum_sq(), other.num(), stats.min(), stats.max()
        )


@pytest.mark.parametrize("protocol", range(pickle.HIGHEST_PROTOCOL + 1))
def test_multi_panel_background_statistics_pickle(protocol):
    multi = MultiPanelImageVolume()
    multi.add(make_image_volume(0, 4, 4, 5))
    multi.add(make_image_volume(0, 4, 6, 3))
    stats = MultiPanelBackgroundStatistics(multi)

    stats2 = pickle.loads(pickle.dumps(stats, protocol=protocol))
    assert len(stats2) == len(stats)
    for i in range(len(stats)):
        assert_statistics_equal(stats.get(i), stats2.get(i))


def test_unpickled_multi_panel_background_statistics_can_be_accumulated():
    def make_multi(frame0, frame1):
        multi = MultiPanelImageVolume()
        multi.add(make_image_volume(frame0, frame1, 4, 5))
        multi.add(make_image_volume(frame0, frame1, 6, 3))
        return MultiPanelBackgroundStatistics(multi)

    stats1, stats2 = make_multi(0, 3), make_multi(3, 6)

    expected = pickle.loads(pickle.dumps(stats1))
    expected += stats2

    result = pickle.loads(pickle.dumps(stats1))
    result += pickle.loads(pickle.dumps(stats2))
    assert len(result) == 2
    for i in range(len(result)):
        assert_statistics_equal(expected.get(i), result.get(i))
