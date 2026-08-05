from __future__ import annotations

import random

import pytest


def test_consistent():
    from dials.array_family import flex
    from dials.model.data import Shoebox

    shoebox = flex.shoebox(10)

    for i in range(10):
        x0 = random.randint(0, 1000)
        y0 = random.randint(0, 1000)
        z0 = random.randint(0, 1000)
        x1 = random.randint(1, 10) + x0
        y1 = random.randint(1, 10) + y0
        z1 = random.randint(1, 10) + z0
        shoebox[i] = Shoebox((x0, x1, y0, y1, z0, z1))

    assert shoebox.is_consistent() == flex.bool(10, False)

    for i in range(10):
        shoebox[i].allocate_data()
        shoebox[i].allocate_background()

    assert shoebox.is_consistent() == flex.bool(10, True)

    for i in [0, 2, 4, 6, 8]:
        shoebox[i].data.resize(flex.grid(10, 10, 10))

    assert shoebox.is_consistent() == flex.bool([False, True] * 5)

    for i in range(10):
        shoebox[i].deallocate()

    assert shoebox.is_consistent() == flex.bool(10, False)


def test_is_bbox_within_image_volume():
    from dials.array_family import flex
    from dials.model.data import Shoebox

    isize = (1000, 1000)
    srange = (0, 100)

    shoebox = flex.shoebox(7)
    shoebox[0] = Shoebox((10, 20, 10, 20, 10, 20))
    shoebox[1] = Shoebox((-10, 20, 10, 20, 10, 20))
    shoebox[2] = Shoebox((10, 20, -10, 20, 10, 20))
    shoebox[3] = Shoebox((10, 20, 10, 20, -10, 20))
    shoebox[4] = Shoebox((10, 1020, 10, 20, 10, 20))
    shoebox[5] = Shoebox((10, 20, 10, 1020, 10, 20))
    shoebox[6] = Shoebox((10, 20, 10, 20, 10, 1020))

    assert shoebox.is_bbox_within_image_volume(isize, srange) == flex.bool(
        [True, False, False, False, False, False, False]
    )


def test_does_bbox_contain_bad_pixels():
    from dials.array_family import flex
    from dials.model.data import Shoebox

    mask = flex.bool(flex.grid(100, 100), True)
    for j in range(100):
        for i in range(40, 60):
            mask[j, i] = False
            mask[i, j] = False

    shoebox = flex.shoebox(1000)
    res = flex.bool(1000)
    for i in range(1000):
        x0 = random.randint(0, 90)
        y0 = random.randint(0, 90)
        z0 = random.randint(0, 90)
        x1 = random.randint(1, 10) + x0
        y1 = random.randint(1, 10) + y0
        z1 = random.randint(1, 10) + z0

        shoebox[i] = Shoebox((x0, x1, y0, y1, z0, z1))

        res2 = False
        if x0 >= 40 and x0 < 60:
            res2 = True
        if x1 > 40 and x1 <= 60:
            res2 = True
        if y0 >= 40 and y0 < 60:
            res2 = True
        if y1 > 40 and y1 <= 60:
            res2 = True

        res[i] = res2

    assert shoebox.does_bbox_contain_bad_pixels(mask) == res


def test_count_mask_values():
    from dials.array_family import flex
    from dials.model.data import Shoebox

    shoebox = flex.shoebox(10)
    num = flex.int(10)
    value = 1 << 2
    for i in range(10):
        x0 = random.randint(0, 90)
        y0 = random.randint(0, 90)
        z0 = random.randint(0, 90)
        x1 = random.randint(1, 10) + x0
        y1 = random.randint(1, 10) + y0
        z1 = random.randint(1, 10) + z0

        shoebox[i] = Shoebox((x0, x1, y0, y1, z0, z1))
        shoebox[i].allocate_data()
        shoebox[i].allocate_background()
        maxnum = len(shoebox[i].mask)
        num[i] = random.randint(1, maxnum)
        indices = random.sample(list(range(maxnum)), num[i])
        for j in indices:
            shoebox[i].mask[j] = value

    assert shoebox.count_mask_values(value) == num


def test_bounding_boxes():
    from dials.array_family import flex
    from dials.model.data import Shoebox

    shoebox = flex.shoebox(10)
    bbox = flex.int6(10)
    for i in range(10):
        x0 = random.randint(0, 90)
        y0 = random.randint(0, 90)
        z0 = random.randint(0, 90)
        x1 = random.randint(1, 10) + x0
        y1 = random.randint(1, 10) + y0
        z1 = random.randint(1, 10) + z0
        bbox[i] = (x0, x1, y0, y1, z0, z1)
        shoebox[i] = Shoebox(bbox[i])

    bbox2 = shoebox.bounding_boxes()
    for i in range(10):
        assert bbox2[i] == bbox[i]


def frame_sliced_shoebox_test_phi():
    """Rotation angles at the centres of images 1 to 13, for a scan of 0.1 radian
    oscillations beginning at zero, alongside the image number of the first"""
    from dials.array_family import flex

    return flex.double([0.1 * (f - 0.5) for f in range(1, 14)]), 1


def frame_sliced_shoebox_test_data():
    """Two shoeboxes with a different number of frames each"""
    from dials.algorithms.shoebox import MaskCode
    from dials.array_family import flex
    from dials.model.data import Shoebox

    shoeboxes = []
    for z0, nz in ((10, 3), (0, 2)):
        shoebox = Shoebox(0, (0, 4, 0, 3, z0, z0 + nz))
        shoebox.allocate_data()
        shoebox.allocate_background()
        for k in range(nz):
            for j in range(3):
                for i in range(4):
                    shoebox.data[k, j, i] = k + 1
                    shoebox.background[k, j, i] = 0.5
                    shoebox.mask[k, j, i] = MaskCode.Valid | MaskCode.Foreground
        shoeboxes.append(shoebox)
    return flex.shoebox(shoeboxes)


def test_frame_sliced_shoebox_from_shoeboxes():
    from dials.array_family import flex

    sliced = flex.frame_sliced_shoebox(
        frame_sliced_shoebox_test_data(), *frame_sliced_shoebox_test_phi()
    )

    assert len(sliced) == 2
    assert list(sliced.num_frames()) == [3, 2]
    assert list(sliced[0].frames) == [11, 12, 13]
    assert list(sliced[1].frames) == [1, 2]
    # Each shoebox takes the rotation angles of its own frames from the scan
    assert list(sliced[0].phi) == pytest.approx([1.05, 1.15, 1.25])
    assert list(sliced[1].phi) == pytest.approx([0.05, 0.15])
    assert list(sliced[0].foreground_sum_raw) == [12.0, 24.0, 36.0]
    assert list(sliced[1].foreground_sum_minus_background) == [6.0, 18.0]

    # Every pixel is a valid foreground pixel, so the summation intensity is the
    # background-subtracted foreground sum, and no frame has a background pixel
    # to contribute an m/n term to the variance
    assert list(sliced[0].summation_intensity) == [6.0, 18.0, 30.0]
    assert list(sliced[0].summation_intensity_variance) == [12.0, 24.0, 36.0]
    assert list(sliced[1].summation_intensity_valid) == [True, True]


def test_frame_sliced_shoebox_data_arrays_round_trip():
    from dials.array_family import flex

    sliced = flex.frame_sliced_shoebox(
        frame_sliced_shoebox_test_data(), *frame_sliced_shoebox_test_phi()
    )
    arrays = sliced.get_frame_sliced_shoebox_data_arrays()

    # The per-frame arrays of every element are concatenated together
    assert list(arrays[0]) == [11, 12, 13, 1, 2]

    rebuilt = flex.frame_sliced_shoebox(sliced.num_frames(), *arrays)
    assert list(rebuilt.num_frames()) == list(sliced.num_frames())
    assert all(a == b for a, b in zip(rebuilt, sliced))


def test_frame_sliced_shoebox_is_picklable():
    import pickle

    from dials.array_family import flex

    sliced = flex.frame_sliced_shoebox(
        frame_sliced_shoebox_test_data(), *frame_sliced_shoebox_test_phi()
    )

    unpickled = pickle.loads(pickle.dumps(sliced, protocol=pickle.HIGHEST_PROTOCOL))
    assert list(unpickled.num_frames()) == [3, 2]
    assert all(a == b for a, b in zip(unpickled, sliced))

    # Default constructed elements have no frames at all
    empty = pickle.loads(
        pickle.dumps(flex.frame_sliced_shoebox(3), protocol=pickle.HIGHEST_PROTOCOL)
    )
    assert list(empty.num_frames()) == [0, 0, 0]
