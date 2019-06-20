from __future__ import absolute_import, division, print_function

import math
import random


def test_getters():
    from dials.algorithms.profile_model.modeller import SingleSampler

    scan_range = (2, 12)
    depth = scan_range[1] - scan_range[0]
    nz = 2
    sampler = SingleSampler(scan_range, nz)
    scan_range = sampler.scan_range()
    grid_size = sampler.grid_size()
    step_size = sampler.step_size()
    size = len(sampler)

    assert scan_range[0] == scan_range[0]
    assert scan_range[1] == scan_range[1]
    assert nz == grid_size
    assert step_size == depth / nz
    assert nz == size


def test_indexing():
    from dials.algorithms.profile_model.modeller import SingleSampler

    width = 1000
    height = 1000
    scan_range = (2, 12)
    nz = 2
    sampler = SingleSampler(scan_range, nz)
    zstep = sampler.step_size()
    zind = [k for k in range(nz)]

    xp = [width * 0.5 for k in range(nz)]
    yp = [height * 0.5 for k in range(nz)]
    zp = [(z + 0.5) * zstep + scan_range[0] for z in zind]

    eps = 1e-10

    for x0, y0, z0, i in zip(xp, yp, zp, range(len(sampler))):
        x1, y1, z1 = sampler.coord(i)
        assert abs(z0 - z1) <= eps


def test_nearest():
    from dials.algorithms.profile_model.modeller import SingleSampler

    scan_range = (2, 12)
    depth = scan_range[1] - scan_range[0]
    nz = 2
    sampler = SingleSampler(scan_range, nz)

    for i in range(1000):
        x = random.uniform(0, 1000)
        y = random.uniform(0, 1000)
        z = random.uniform(*scan_range)
        k = int(math.floor((z - scan_range[0]) / (depth / nz)))
        if k >= nz:
            k = nz - 1
        index0 = k
        index1 = sampler.nearest(0, (x, y, z))
        assert index0 == index1


def test_nearest_n():
    from dials.algorithms.profile_model.modeller import SingleSampler

    scan_range = (2, 12)
    depth = scan_range[1] - scan_range[0]
    nz = 2
    sampler = SingleSampler(scan_range, nz)

    for i in range(1000):
        x = random.uniform(0, 1000)
        y = random.uniform(0, 1000)
        z = random.uniform(*scan_range)
        k = int(math.floor((z - scan_range[0]) * nz / depth))
        if k >= nz:
            k = nz - 1
        index0 = k
        index1 = sampler.nearest_n(0, (x, y, z))
        assert len(set(index1)) == len(index1)
        assert index0 == index1[-1]
        for ind in index1:
            assert abs(ind - k) <= 1


def test_weights():
    from dials.algorithms.profile_model.modeller import SingleSampler
    from scitbx import matrix

    scan_range = (2, 12)
    nz = 2
    sampler = SingleSampler(scan_range, nz)

    # Check the weight at the coord in 1.0
    eps = 1e-7
    for i in range(len(sampler)):
        coord = sampler.coord(i)
        weight = sampler.weight(i, 0, coord)
        assert abs(weight - 1.0) < eps

    # Ensure we get the expected weight at the next grid point at half way
    # between grid points
    expected = math.exp(-4.0 * math.log(2.0))
    for k in range(nz):
        coord1 = matrix.col(sampler.coord(k))
        if k > 0:
            coord = matrix.col(sampler.coord(k - 1))
            weight = sampler.weight(k, 0, coord)
            assert abs(weight - expected) < eps
            weight = sampler.weight(k, 0, (coord + coord1) / 2.0)
            assert abs(weight - 0.5) < eps
        if k < nz - 1:
            coord = matrix.col(sampler.coord(k + 1))
            weight = sampler.weight(k, 0, coord)
            assert abs(weight - expected) < eps
            weight = sampler.weight(k, 0, (coord + coord1) / 2.0)
            assert abs(weight - 0.5) < eps


def test_self_consistent():
    from dials.algorithms.profile_model.modeller import SingleSampler

    scan_range = (2, 12)
    nz = 2
    sampler = SingleSampler(scan_range, nz)

    for i in range(len(sampler)):
        coord = sampler.coord(i)
        index = sampler.nearest(0, coord)
        assert index == i


def test_pickle():
    from dials.algorithms.profile_model.modeller import SingleSampler
    import six.moves.cPickle as pickle

    scan_range = (2, 12)
    nz = 2
    sampler = SingleSampler(scan_range, nz)
    sampler2 = pickle.loads(pickle.dumps(sampler))
    assert sampler.grid_size() == sampler2.grid_size()
