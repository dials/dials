from __future__ import annotations

import random


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


def test_peak_coordinates_unique_maximum():
    from dials.array_family import flex
    from dials.model.data import Shoebox

    # A shoebox with a single, unique brightest pixel: the peak must be that
    # pixel's global coordinate (bbox origin + index + 0.5). This covers the
    # cheap, untied path which must behave exactly as before.
    bbox = (10, 14, 20, 24, 30, 34)
    shoebox = flex.shoebox(1)
    shoebox[0] = Shoebox(bbox)
    shoebox[0].allocate_data()
    data = shoebox[0].data  # shape (z, y, x)
    data[1, 2, 3] = 100

    peak = shoebox.peak_coordinates()
    assert peak[0] == (10 + 3 + 0.5, 20 + 2 + 0.5, 30 + 1 + 0.5)


def test_peak_coordinates_ties_resolved_by_centroid():
    from dials.array_family import flex
    from dials.model.data import Shoebox

    # Two pixels share the maximum value. Place the intensity mass (and hence the
    # centroid) near one of them; the tie must be broken in favour of the pixel
    # nearest the centroid rather than the z, y, x-first pixel. On the previous
    # (max_index only) implementation this returned the far corner instead.
    # See https://github.com/dials/dials/issues/3014.
    bbox = (0, 6, 0, 6, 0, 1)  # single frame, 6x6 in y, x
    shoebox = flex.shoebox(1)
    shoebox[0] = Shoebox(bbox)
    shoebox[0].allocate_data()
    data = shoebox[0].data  # shape (z=1, y=6, x=6)

    # Build a blob of intensity in the high-y, high-x corner so the centroid sits
    # there. Two equal maxima: one at the near corner (in the blob) and one at
    # the far (low-y, low-x) corner which is first in z, y, x order.
    for y in range(3, 6):
        for x in range(3, 6):
            data[0, y, x] = 5
    data[0, 0, 0] = 10  # far, z/y/x-first maximum
    data[0, 5, 5] = 10  # near the centroid

    peak = shoebox.peak_coordinates()
    # Expect the near-centroid pixel (x=5, y=5), not the (0, 0) corner.
    assert peak[0] == (5 + 0.5, 5 + 0.5, 0 + 0.5)


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
