from __future__ import absolute_import, division, print_function

import math
import random

from dials.model.data import Shoebox
from scitbx import matrix
import six.moves.cPickle as pickle


def random_shoeboxes(num, mask=False):
    for i in range(num):
        x0 = random.randint(0, 100)
        y0 = random.randint(0, 100)
        z0 = random.randint(0, 100)
        x1 = random.randint(x0 + 5, x0 + 20)
        y1 = random.randint(y0 + 5, y0 + 20)
        z1 = random.randint(z0 + 5, z0 + 20)
        bbox = (x0, x1, y0, y1, z0, z1)
        xc0 = (x1 + x0) / 2.0
        yc0 = (y1 + y0) / 2.0
        zc0 = (z1 + z0) / 2.0
        xc = random.uniform(xc0 - 1, xc0 + 1)
        yc = random.uniform(yc0 - 1, yc0 + 1)
        zc = random.uniform(zc0 - 1, zc0 + 1)
        centre = (xc, yc, zc)
        intensity = random.randint(10, 10000)
        shoebox = generate_shoebox(bbox, centre, intensity, mask=mask)
        yield (shoebox, (centre, intensity))


def generate_shoebox(bbox, centre, intensity, mask=False):
    from dials.algorithms.shoebox import MaskCode

    shoebox = Shoebox()
    shoebox.bbox = bbox
    shoebox.allocate()
    for i in range(len(shoebox.mask)):
        shoebox.mask[i] = MaskCode.Valid | MaskCode.Foreground
    shoebox.data = gaussian(
        shoebox.size(),
        1.0,
        [c - o for c, o in zip(centre[::-1], shoebox.offset())],
        [s / 8.0 for s in shoebox.size()],
    )
    if mask:
        shoebox.mask = create_mask(
            shoebox.size(),
            [c - o for c, o in zip(centre[::-1], shoebox.offset())],
            MaskCode.Valid | MaskCode.Foreground,
        )
    tot = 0
    mask_code = MaskCode.Valid | MaskCode.Foreground
    for i in range(len(shoebox.data)):
        if shoebox.mask[i] & mask_code == mask_code:
            tot += shoebox.data[i]
    if tot > 0:
        shoebox.data *= intensity / tot
    return shoebox


def create_mask(size, x0, value):
    from scitbx.array_family import flex

    mask = flex.int(flex.grid(size), 0)
    rad = min(s - c for s, c in zip(size, x0))
    for k in range(size[0]):
        for j in range(size[1]):
            for i in range(size[2]):
                d = math.sqrt((j - x0[1]) ** 2 + (i - x0[2]) ** 2)
                if d < rad:
                    mask[k, j, i] = value
    return mask


def evaluate_gaussian(x, a, x0, sx):
    assert len(x) == len(x0)
    assert len(x) == len(sx)

    g = 0.0
    for xi, x0i, sxi in zip(x, x0, sx):
        g += (xi - x0i) ** 2 / (2.0 * sxi ** 2)

    return a * math.exp(-g)


def gaussian(size, a, x0, sx):
    from dials.array_family import flex

    result = flex.real(flex.grid(size))

    index = [0] * len(size)
    while True:
        result[index] = evaluate_gaussian(index, a, x0, sx)
        for j in range(len(size)):
            index[j] += 1
            if index[j] < size[j]:
                break
            index[j] = 0
            if j == len(size) - 1:
                return result


def test_allocate():
    for i in range(10):
        x0 = random.randint(0, 1000)
        y0 = random.randint(0, 1000)
        z0 = random.randint(0, 1000)
        x1 = random.randint(1, 10) + x0
        y1 = random.randint(1, 10) + y0
        z1 = random.randint(1, 10) + z0

        shoebox = Shoebox((x0, x1, y0, y1, z0, z1))
        shoebox.allocate()
        assert shoebox.data.all() == (z1 - z0, y1 - y0, x1 - x0)
        assert shoebox.mask.all() == (z1 - z0, y1 - y0, x1 - x0)
        shoebox.deallocate()
        assert shoebox.data.all() == (0, 0, 0)
        assert shoebox.mask.all() == (0, 0, 0)


def test_offset():
    for i in range(10):
        x0 = random.randint(0, 1000)
        y0 = random.randint(0, 1000)
        z0 = random.randint(0, 1000)
        x1 = random.randint(1, 10) + x0
        y1 = random.randint(1, 10) + y0
        z1 = random.randint(1, 10) + z0

        shoebox = Shoebox((x0, x1, y0, y1, z0, z1))
        assert shoebox.xoffset() == x0
        assert shoebox.yoffset() == y0
        assert shoebox.zoffset() == z0
        assert shoebox.offset() == (z0, y0, x0)


def test_size():
    for i in range(10):
        x0 = random.randint(0, 1000)
        y0 = random.randint(0, 1000)
        z0 = random.randint(0, 1000)
        x1 = random.randint(1, 10) + x0
        y1 = random.randint(1, 10) + y0
        z1 = random.randint(1, 10) + z0

        shoebox = Shoebox((x0, x1, y0, y1, z0, z1))
        assert shoebox.xsize() == x1 - x0
        assert shoebox.ysize() == y1 - y0
        assert shoebox.zsize() == z1 - z0
        assert shoebox.size() == (z1 - z0, y1 - y0, x1 - x0)


def test_consistent():
    from dials.array_family import flex

    for i in range(1000):
        x0 = random.randint(0, 1000)
        y0 = random.randint(0, 1000)
        z0 = random.randint(0, 1000)
        x1 = random.randint(1, 10) + x0
        y1 = random.randint(1, 10) + y0
        z1 = random.randint(1, 10) + z0
        try:
            shoebox = Shoebox((x0, x1, y0, y1, z0, z1))
            assert not shoebox.is_consistent()
            shoebox.allocate()
            assert shoebox.is_consistent()
            shoebox.data = flex.real(flex.grid(20, 20, 20))
            assert not shoebox.is_consistent()
            shoebox.deallocate()
            assert not shoebox.is_consistent()
        except Exception:
            print(x0, y0, z0, x1, y1, z1)
            raise


def test_is_bbox_within_image_volume():
    isize = (1000, 1000)
    srange = (0, 100)

    shoebox = Shoebox((10, 20, 10, 20, 10, 20))
    assert shoebox.is_bbox_within_image_volume(isize, srange)
    shoebox = Shoebox((-10, 20, 10, 20, 10, 20))
    assert not shoebox.is_bbox_within_image_volume(isize, srange)
    shoebox = Shoebox((10, 20, -10, 20, 10, 20))
    assert not shoebox.is_bbox_within_image_volume(isize, srange)
    shoebox = Shoebox((10, 20, 10, 20, -10, 20))
    assert not shoebox.is_bbox_within_image_volume(isize, srange)
    shoebox = Shoebox((10, 1020, 10, 20, 10, 20))
    assert not shoebox.is_bbox_within_image_volume(isize, srange)
    shoebox = Shoebox((10, 20, 10, 1020, 10, 20))
    assert not shoebox.is_bbox_within_image_volume(isize, srange)
    shoebox = Shoebox((10, 20, 10, 20, 10, 1020))
    assert not shoebox.is_bbox_within_image_volume(isize, srange)


def test_does_bbox_contain_bad_pixels():
    from scitbx.array_family import flex

    mask = flex.bool(flex.grid(100, 100), True)
    for j in range(100):
        for i in range(40, 60):
            mask[j, i] = False
            mask[i, j] = False

    for i in range(1000):
        x0 = random.randint(0, 90)
        y0 = random.randint(0, 90)
        z0 = random.randint(0, 90)
        x1 = random.randint(1, 10) + x0
        y1 = random.randint(1, 10) + y0
        z1 = random.randint(1, 10) + z0

        shoebox = Shoebox((x0, x1, y0, y1, z0, z1))

        res1 = shoebox.does_bbox_contain_bad_pixels(mask)
        res2 = False
        if x0 >= 40 and x0 < 60:
            res2 = True
        if x1 > 40 and x1 <= 60:
            res2 = True
        if y0 >= 40 and y0 < 60:
            res2 = True
        if y1 > 40 and y1 <= 60:
            res2 = True

        assert res1 == res2


def test_count_mask_values():
    for i in range(10):
        x0 = random.randint(0, 90)
        y0 = random.randint(0, 90)
        z0 = random.randint(0, 90)
        x1 = random.randint(1, 10) + x0
        y1 = random.randint(1, 10) + y0
        z1 = random.randint(1, 10) + z0

        shoebox = Shoebox((x0, x1, y0, y1, z0, z1))
        shoebox.allocate()
        maxnum = len(shoebox.mask)
        num = random.randint(1, maxnum)
        indices = random.sample(list(range(maxnum)), num)
        value = 1 << 2
        for i in indices:
            shoebox.mask[i] = value

        assert shoebox.count_mask_values(value) == num


def test_centroid_all():
    for shoebox, (XC, I) in random_shoeboxes(10):
        centroid = shoebox.centroid_all()
        assert shoebox.is_consistent()
        assert abs(matrix.col(centroid.px.position) - matrix.col(XC)) < 1.0


def test_centroid_masked():
    for shoebox, (XC, I) in random_shoeboxes(10):
        centroid = shoebox.centroid_masked(1 << 0)
        assert shoebox.is_consistent()
        assert abs(matrix.col(centroid.px.position) - matrix.col(XC)) < 1.0


def test_summed_intensity():
    for shoebox, (XC, I) in random_shoeboxes(10):
        intensity = shoebox.summed_intensity()
        assert shoebox.is_consistent()
        assert abs(intensity.observed.value - I) < 1e-1


def test_flatten():
    from dials.array_family import flex
    from dials.algorithms.shoebox import MaskCode

    for shoebox, (XC, I) in random_shoeboxes(10, mask=True):
        assert not shoebox.flat
        zs = shoebox.zsize()
        ys = shoebox.ysize()
        xs = shoebox.xsize()
        expected_data = flex.real(flex.grid(1, ys, xs), 0)
        expected_mask = flex.int(flex.grid(1, ys, xs), 0)
        for k in range(zs):
            for j in range(ys):
                for i in range(xs):
                    expected_data[0, j, i] += shoebox.data[k, j, i]
                    expected_mask[0, j, i] |= shoebox.mask[k, j, i]
                    if not (expected_mask[0, j, i] & MaskCode.Valid) or not (
                        shoebox.mask[k, j, i] & MaskCode.Valid
                    ):
                        expected_mask[0, j, i] &= ~MaskCode.Valid
        shoebox.flatten()
        diff = expected_data.as_double() - shoebox.data.as_double()
        max_diff = flex.max(flex.abs(diff))
        assert max_diff < 1e-7
        assert expected_mask.all_eq(shoebox.mask)
        assert shoebox.flat
        assert shoebox.is_consistent()


def test_all_foreground_valid():
    from dials.test.model.data.all_foreground_valid_data import data

    shoeboxes = pickle.loads(data)
    for i, shoebox in enumerate(shoeboxes):
        if i < 4:
            assert not shoebox.all_foreground_valid()
        else:
            assert shoebox.all_foreground_valid()
