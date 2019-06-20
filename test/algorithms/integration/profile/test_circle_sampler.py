from __future__ import absolute_import, division, print_function

import math
import random


def test_getters():
    from dials.algorithms.profile_model.modeller import CircleSampler

    width = 1000
    height = 1000
    scan_range = (2, 12)
    nz = 2
    sampler = CircleSampler((width, height), scan_range, nz)
    image_size = sampler.image_size()
    image_centre = sampler.image_centre()
    r0 = sampler.r0()
    r1 = sampler.r1()
    r2 = sampler.r2()
    size = len(sampler)

    assert width == image_size[0]
    assert height == image_size[1]
    assert width // 2 == image_centre[0]
    assert height // 2 == image_centre[1]
    assert r0 == min([width // 2, height // 2])
    assert r1 == r0 / 3.0
    assert r2 == r1 * math.sqrt(5.0)
    assert 9 * nz == size


def test_detector_area():
    from dials.algorithms.profile_model.modeller import CircleSampler
    from scitbx.array_family import flex

    width = 1000
    height = 1000
    scan_range = (2, 12)
    nz = 2
    sampler = CircleSampler((width, height), scan_range, nz)
    im = flex.int(flex.grid(height, width))
    for j in range(height):
        for i in range(width):
            im[j, i] = sampler.nearest(0, (i, j, 0))

    assert im[height // 2, width // 2] == 0
    assert im[height // 2, width - 1] == 1
    assert im[height - 1, width - 1] == 2
    assert im[height - 1, width // 2] == 3
    assert im[height - 1, 0] == 4
    assert im[height // 2, 0] == 5
    assert im[0, 0] == 6
    assert im[0, width // 2] == 7
    assert im[0, width - 1] == 8


#        from matplotlib import pylab
#        pylab.imshow(im.as_numpy_array())
#        pylab.show()


def test_indexing():
    from dials.algorithms.profile_model.modeller import CircleSampler

    width = 1000
    height = 1000
    scan_range = (2, 12)
    depth = scan_range[1] - scan_range[0]
    nz = 2
    sampler = CircleSampler((width, height), scan_range, nz)
    zstep = depth / nz
    eps = 1e-10

    dt = 2 * math.pi / 8
    r2 = sampler.r2()

    xp = [width / 2.0] + [width / 2.0 + r2 * math.cos(i * dt) for i in range(8)]
    yp = [height / 2.0] + [height / 2.0 + r2 * math.sin(i * dt) for i in range(8)]
    zind = [[k] * 9 for k in range(nz)]
    zind = [i for j in zind for i in j]
    zp = [(z + 0.5) * zstep + scan_range[0] for z in zind]

    for x0, y0, z0, i in zip(xp, yp, zp, range(len(sampler))):
        x1, y1, z1 = sampler.coord(i)
        assert abs(x0 - x1) <= eps
        assert abs(y0 - y1) <= eps
        assert abs(z0 - z1) <= eps


def test_nearest():
    from dials.algorithms.profile_model.modeller import CircleSampler

    width = 1000
    height = 1000
    scan_range = (2, 12)
    depth = scan_range[1] - scan_range[0]
    nz = 2
    sampler = CircleSampler((width, height), scan_range, nz)
    xc, yc = sampler.image_centre()
    r1 = sampler.r1()

    for i in range(1000):
        x = random.uniform(0, 1000)
        y = random.uniform(0, 1000)
        z = random.uniform(scan_range[0], scan_range[1])

        r = math.sqrt((x - xc) ** 2 + (y - yc) ** 2)
        if r < r1:
            index00 = 0
        else:
            t = math.atan2(y - yc, x - xc)
            ai = int(math.floor(t * 8 / (2 * math.pi) + 0.5)) % 8
            index00 = ai + 1

        index01 = int((z - scan_range[0]) * nz / depth)
        if index01 < 0:
            index01 = 0
        if index01 >= 2:
            index01 = 1

        index0 = index00 + index01 * 9
        index1 = sampler.nearest(0, (x, y, z))
        assert index0 == index1


def test_nearest_n():
    from dials.algorithms.profile_model.modeller import CircleSampler

    width = 1000
    height = 1000
    scan_range = (2, 12)
    depth = scan_range[1] - scan_range[0]
    nz = 2
    sampler = CircleSampler((width, height), scan_range, nz)
    xc, yc = sampler.image_centre()
    r1 = sampler.r1()

    for i in range(1000):
        x = random.uniform(0, 1000)
        y = random.uniform(0, 1000)
        z = random.uniform(scan_range[0], scan_range[1])

        r = math.sqrt((x - xc) ** 2 + (y - yc) ** 2)
        if r < r1:
            index00 = 0
        else:
            t = math.atan2(y - yc, x - xc)
            ai = int(math.floor(t * 8 / (2 * math.pi) + 0.5)) % 8
            index00 = ai + 1

        index01 = int((z - scan_range[0]) * nz / depth)
        if index01 < 0:
            index01 = 0
        if index01 >= 2:
            index01 = 1

        index0 = index00 + index01 * 9
        index1 = sampler.nearest_n(0, (x, y, z))
        assert index0 == index1[0]
        if index0 % 9 == 0:
            assert len(index1) == 9
            assert all(idx == index0 + i for i, idx in enumerate(index1))
        else:
            assert len(index1) == 4
            assert index1[1] == (index0 // 9) * 9
            if (index0 % 9) == 1:
                assert index1[2] == index0 + 1
                assert index1[3] == index0 + 7
            elif (index0 % 9) == 8:
                assert index1[2] == index0 - 7
                assert index1[3] == index0 - 1
            else:
                assert index1[2] == index0 + 1
                assert index1[3] == index0 - 1


def test_weights():
    from dials.algorithms.profile_model.modeller import CircleSampler
    from scitbx import matrix

    width = 1000
    height = 1000
    scan_range = (2, 12)
    nz = 1
    sampler = CircleSampler((width, height), scan_range, nz)

    # Check the weight at the coord in 1.0
    eps = 1e-7
    for i in range(len(sampler)):
        coord = sampler.coord(i)
        weight = sampler.weight(i, 0, coord)
        assert abs(weight - 1.0) < eps

    r1 = sampler.r1()
    r2 = sampler.r2()
    r = r2 / (2.0 * r1)
    expected = math.exp(-4.0 * r * r * math.log(2.0))
    for i in range(1, 9):
        coord = sampler.coord(i)
        weight = sampler.weight(0, 0, coord)
        assert abs(weight - expected) < eps

    r = r2 / (2.0 * (r2 - r1))
    expected = math.exp(-4.0 * r * r * math.log(2.0))
    for i in range(1, 9):
        coord = sampler.coord(0)
        weight = sampler.weight(i, 0, coord)
        assert abs(weight - expected) < eps

    for i in range(1, 9):
        coord1 = matrix.col(sampler.coord(0))
        coord2 = matrix.col(sampler.coord(i))
        coord = coord1 + r1 * (coord2 - coord1) / r2
        weight = sampler.weight(0, 0, coord)
        assert abs(weight - 0.5) < eps
        weight = sampler.weight(i, 0, coord)
        assert abs(weight - 0.5) < eps


def test_self_consistent():
    from dials.algorithms.profile_model.modeller import CircleSampler

    width = 1000
    height = 1000
    scan_range = (2, 12)
    nz = 2
    sampler = CircleSampler((width, height), scan_range, nz)

    for i in range(len(sampler)):
        coord = sampler.coord(i)
        index = sampler.nearest(0, coord)
        assert index == i


def test_z_index():
    from dials.algorithms.profile_model.modeller import CircleSampler

    width = 1000
    height = 1000
    scan_range = (2, 12)
    nz = 2
    sampler = CircleSampler((width, height), scan_range, nz)
    assert (sampler.nearest(0, (500, 500, 2.0)) / 9) == 0
    assert (sampler.nearest(0, (500, 500, 3.0)) / 9) == 0
    assert (sampler.nearest(0, (500, 500, 4.0)) / 9) == 0
    assert (sampler.nearest(0, (500, 500, 5.0)) / 9) == 0
    assert (sampler.nearest(0, (500, 500, 6.0)) / 9) == 0
    assert (sampler.nearest(0, (500, 500, 6.5)) / 9) == 0
    assert (sampler.nearest(0, (500, 500, 7.0)) / 9) == 1
    assert (sampler.nearest(0, (500, 500, 7.5)) / 9) == 1
    assert (sampler.nearest(0, (500, 500, 8.0)) / 9) == 1
    assert (sampler.nearest(0, (500, 500, 9.0)) / 9) == 1
    assert (sampler.nearest(0, (500, 500, 10.0)) / 9) == 1
    assert (sampler.nearest(0, (500, 500, 11.0)) / 9) == 1
    assert (sampler.nearest(0, (500, 500, 12.0)) / 9) == 1


def test_pickle():
    from dials.algorithms.profile_model.modeller import CircleSampler
    import six.moves.cPickle as pickle

    width = 1000
    height = 1000
    scan_range = (2, 12)
    nz = 2
    sampler = CircleSampler((width, height), scan_range, nz)

    sampler2 = pickle.loads(pickle.dumps(sampler))

    assert sampler.image_size() == sampler2.image_size()
    assert sampler.num_z() == sampler2.num_z()
