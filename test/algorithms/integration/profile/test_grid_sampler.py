from __future__ import absolute_import, division, print_function

import math
import random


def test_getters():
    from dials.algorithms.profile_model.modeller import GridSampler

    width = 1000
    height = 1000
    scan_range = (2, 12)
    depth = scan_range[1] - scan_range[0]
    nx = 10
    ny = 10
    nz = 2
    sampler = GridSampler((width, height), scan_range, (nx, ny, nz))
    image_size = sampler.image_size()
    scan_range = sampler.scan_range()
    grid_size = sampler.grid_size()
    step_size = sampler.step_size()
    size = len(sampler)

    assert width == image_size[0]
    assert height == image_size[1]
    assert scan_range[0] == scan_range[0]
    assert scan_range[1] == scan_range[1]
    assert nx == grid_size[0]
    assert ny == grid_size[1]
    assert nz == grid_size[2]
    assert step_size[0] == width / nx
    assert step_size[1] == height / ny
    assert step_size[2] == depth / nz
    assert nx * ny * nz == size


def test_indexing():
    from dials.algorithms.profile_model.modeller import GridSampler

    width = 1000
    height = 1000
    scan_range = (2, 12)
    nx = 10
    ny = 10
    nz = 2
    sampler = GridSampler((width, height), scan_range, (nx, ny, nz))
    xstep, ystep, zstep = sampler.step_size()
    xind = [[i for i in range(nx)]] * ny * nz
    yind = [[j] * nx for j in range(ny)] * nz
    zind = [[k] * nx * ny for k in range(nz)]
    xind = [i for j in xind for i in j]
    yind = [i for j in yind for i in j]
    zind = [i for j in zind for i in j]

    xp = [(x + 0.5) * xstep for x in xind]
    yp = [(y + 0.5) * ystep for y in yind]
    zp = [(z + 0.5) * zstep + scan_range[0] for z in zind]

    eps = 1e-10

    for x0, y0, z0, i in zip(xp, yp, zp, range(len(sampler))):
        x1, y1, z1 = sampler.coord(i)
        assert abs(x0 - x1) <= eps
        assert abs(y0 - y1) <= eps
        assert abs(z0 - z1) <= eps


def test_nearest():
    from dials.algorithms.profile_model.modeller import GridSampler

    width = 1000
    height = 1000
    scan_range = (2, 12)
    depth = scan_range[1] - scan_range[0]
    nx = 10
    ny = 10
    nz = 2
    sampler = GridSampler((width, height), scan_range, (nx, ny, nz))

    for i in range(1000):
        x = random.uniform(0, 1000)
        y = random.uniform(0, 1000)
        z = random.uniform(*scan_range)
        i = int(math.floor(x / (width / nx)))
        j = int(math.floor(y / (height / ny)))
        k = int(math.floor((z - scan_range[0]) / (depth / nz)))
        if i >= nx:
            i = nx - 1
        if j >= ny:
            j = ny - 1
        if k >= nz:
            k = nz - 1
        index0 = i + j * nx + k * nx * ny
        index1 = sampler.nearest(0, (x, y, z))
        assert index0 == index1


def test_nearest_n():
    from dials.algorithms.profile_model.modeller import GridSampler

    width = 1000
    height = 1000
    scan_range = (2, 12)
    depth = scan_range[1] - scan_range[0]
    nx = 10
    ny = 10
    nz = 2
    sampler = GridSampler((width, height), scan_range, (nx, ny, nz))

    for i in range(1000):
        x = random.uniform(0, 1000)
        y = random.uniform(0, 1000)
        z = random.uniform(*scan_range)
        i = int(math.floor(x * nx / width))
        j = int(math.floor(y * ny / height))
        k = int(math.floor((z - scan_range[0]) * nz / depth))
        if i >= nx:
            i = nx - 1
        if j >= ny:
            j = ny - 1
        if k >= nz:
            k = nz - 1
        index0 = i + j * nx + k * nx * ny
        index1 = sampler.nearest_n(0, (x, y, z))
        assert len(set(index1)) == len(index1)
        assert index0 == index1[-1]
        for ind in index1:
            ii = ind % nx
            jk = ind // nx
            jj = jk % ny
            kk = jk // ny
            assert abs(ii - i) <= 1
            assert abs(jj - j) <= 1
            assert abs(kk - k) <= 1


def test_weights():
    from dials.algorithms.profile_model.modeller import GridSampler
    from scitbx import matrix

    width = 1000
    height = 1000
    scan_range = (2, 12)
    nx = 10
    ny = 10
    nz = 2
    sampler = GridSampler((width, height), scan_range, (nx, ny, nz))

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
        for j in range(ny):
            for i in range(nx):
                l1 = (i + 0) + ((j + 0) + (k + 0) * ny) * nx
                l2 = (i + 1) + ((j + 0) + (k + 0) * ny) * nx
                l3 = (i - 1) + ((j + 0) + (k + 0) * ny) * nx
                l4 = (i + 0) + ((j + 1) + (k + 0) * ny) * nx
                l5 = (i + 0) + ((j - 1) + (k + 0) * ny) * nx
                l6 = (i + 0) + ((j + 0) + (k + 1) * ny) * nx
                l7 = (i + 0) + ((j + 0) + (k - 1) * ny) * nx
                coord1 = matrix.col(sampler.coord(l1))
                if i < nx - 1:
                    coord = matrix.col(sampler.coord(l2))
                    weight = sampler.weight(l1, 0, coord)
                    assert abs(weight - expected) < eps
                    weight = sampler.weight(l1, 0, (coord + coord1) / 2.0)
                    assert abs(weight - 0.5) < eps
                if i > 0:
                    coord = matrix.col(sampler.coord(l3))
                    weight = sampler.weight(l1, 0, coord)
                    assert abs(weight - expected) < eps
                    weight = sampler.weight(l1, 0, (coord1 + coord) / 2.0)
                    assert abs(weight - 0.5) < eps
                if j < ny - 1:
                    coord = matrix.col(sampler.coord(l4))
                    weight = sampler.weight(l1, 0, coord)
                    assert abs(weight - expected) < eps
                    weight = sampler.weight(l1, 0, (coord + coord1) / 2.0)
                    assert abs(weight - 0.5) < eps
                if j > 0:
                    coord = matrix.col(sampler.coord(l5))
                    weight = sampler.weight(l1, 0, coord)
                    assert abs(weight - expected) < eps
                    weight = sampler.weight(l1, 0, (coord1 + coord) / 2.0)
                    assert abs(weight - 0.5) < eps
                if k < nz - 1:
                    coord = matrix.col(sampler.coord(l6))
                    weight = sampler.weight(l1, 0, coord)
                    assert abs(weight - expected) < eps
                    weight = sampler.weight(l1, 0, (coord + coord1) / 2.0)
                    assert abs(weight - 0.5) < eps
                if k > 0:
                    coord = matrix.col(sampler.coord(l7))
                    weight = sampler.weight(l1, 0, coord)
                    assert abs(weight - expected) < eps
                    weight = sampler.weight(l1, 0, (coord1 + coord) / 2.0)
                    assert abs(weight - 0.5) < eps


def test_self_consistent():
    from dials.algorithms.profile_model.modeller import GridSampler

    width = 1000
    height = 1000
    scan_range = (2, 12)
    nx = 10
    ny = 10
    nz = 2
    sampler = GridSampler((width, height), scan_range, (nx, ny, nz))

    for i in range(len(sampler)):
        coord = sampler.coord(i)
        index = sampler.nearest(0, coord)
        assert index == i


def test_pickle():
    from dials.algorithms.profile_model.modeller import GridSampler
    import six.moves.cPickle as pickle

    width = 1000
    height = 1000
    scan_range = (2, 12)
    nx = 10
    ny = 10
    nz = 2
    sampler = GridSampler((width, height), scan_range, (nx, ny, nz))

    sampler2 = pickle.loads(pickle.dumps(sampler))

    assert sampler.image_size() == sampler2.image_size()
    assert sampler.grid_size() == sampler2.grid_size()
