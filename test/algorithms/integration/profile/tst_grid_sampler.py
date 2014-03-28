from __future__ import division

class Test(object):

  def __init__(self):
    pass

  def run(self):
    self.tst_getters()
    self.tst_indexing()
    self.tst_nearest()
    self.tst_nearest_n()
    self.tst_weights()
    self.tst_self_consistent()
    self.tst_pickle()

  def tst_getters(self):
    from dials.algorithms.integration.profile import GridSampler
    width = 1000
    height = 1000
    depth = 10
    nx = 10
    ny = 10
    nz = 2
    sampler = GridSampler((width, height, depth), (nx, ny, nz))
    volume_size = sampler.volume_size()
    grid_size = sampler.grid_size()
    step_size = sampler.step_size()
    size = len(sampler)

    assert(width == volume_size[0])
    assert(height == volume_size[1])
    assert(depth == volume_size[2])
    assert(nx == grid_size[0])
    assert(ny == grid_size[1])
    assert(nz == grid_size[2])
    assert(step_size[0] == width / nx)
    assert(step_size[1] == height / ny)
    assert(step_size[2] == depth / nz)
    assert(nx * ny * nz == size)
    print 'OK'

  def tst_indexing(self):
    from dials.algorithms.integration.profile import GridSampler
    width = 1000
    height = 1000
    depth = 10
    nx = 10
    ny = 10
    nz = 2
    sampler = GridSampler((width, height, depth), (nx, ny, nz))
    xstep, ystep, zstep = sampler.step_size()
    xind = [[i for i in range(nx)]] * ny * nz
    yind = [[j] * nx for j in range(ny)] * nz
    zind = [[k] * nx * ny for k in range(nz)]
    xind = [i for j in xind for i in j]
    yind = [i for j in yind for i in j]
    zind = [i for j in zind for i in j]

    xp = [(x + 0.5) * xstep for x in xind]
    yp = [(y + 0.5) * ystep for y in yind]
    zp = [(z + 0.5) * zstep for z in zind]

    eps = 1e-10

    for x0, y0, z0, (x1, y1, z1) in zip(xp, yp, zp, sampler):
      assert(abs(x0 - x1) <= eps)
      assert(abs(y0 - y1) <= eps)
      assert(abs(z0 - z1) <= eps)

    print 'OK'


  def tst_nearest(self):
    from random import randint
    from dials.algorithms.integration.profile import GridSampler
    width = 1000
    height = 1000
    depth = 10
    nx = 10
    ny = 10
    nz = 2
    sampler = GridSampler((width, height, depth), (nx, ny, nz))

    for i in range(1000):
      x = randint(0, 1000)
      y = randint(0, 1000)
      z = randint(0, 10)
      i = int((x+0.5) * nx // 1000)
      j = int((y+0.5) * ny // 1000)
      k = int((z+0.5) * nz // 10)
      if i >= nx:
        i = nx - 1
      if j >= ny:
        j = ny - 1
      if k >= nz:
        k = nz - 1
      index0 = i + j * nx + k * nx * ny
      index1 = sampler.nearest((x, y, z))
      assert(index0 == index1)

    print 'OK'

  def tst_nearest_n(self):
    from random import randint
    from dials.algorithms.integration.profile import GridSampler
    width = 1000
    height = 1000
    depth = 10
    nx = 10
    ny = 10
    nz = 2
    sampler = GridSampler((width, height, depth), (nx, ny, nz))

    for i in range(1000):
      x = randint(0, 1000)
      y = randint(0, 1000)
      z = randint(0, 10)
      i = int((x+0.5) * nx // 1000)
      j = int((y+0.5) * ny // 1000)
      k = int((z+0.5) * nz // 10)
      if i >= nx:
        i = nx - 1
      if j >= ny:
        j = ny - 1
      if k >= nz:
        k = nz - 1
      index0 = i + j * nx + k * nx * ny
      index1 = sampler.nearest_n((x, y, z))
      assert(index0 == index1[0])
      c = 1
      if i > 0:
        assert(index1[c] == index0 - 1)
        c += 1
      if j > 0:
        assert(index1[c] == index0 - nx)
        c += 1
      if k > 0:
        assert(index1[c] == index0 - nx * ny)
        c += 1
      if i < nx-1:
        assert(index1[c] == index0 + 1)
        c += 1
      if j < ny-1:
        assert(index1[c] == index0 + nx)
        c += 1
      if k < nz-1:
        assert(index1[c] == index0 + nx * ny)
        c += 1
      assert(c == len(index1))

    print 'OK'

  def tst_weights(self):
    from dials.algorithms.integration.profile import GridSampler
    from scitbx import matrix
    from math import log, exp
    width = 1000
    height = 1000
    depth = 10
    nx = 10
    ny = 10
    nz = 2
    sampler = GridSampler((width, height, depth), (nx, ny, nz))

    # Check the weight at the coord in 1.0
    eps = 1e-7
    for i in range(len(sampler)):
      coord = sampler[i]
      weight = sampler.weight(i, coord)
      assert(abs(weight - 1.0) < eps)

    # Ensure we get the expected weight at the next grid point at half way
    # between grid points
    expected = exp(-4.0 * log(2.0))
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
          coord1 = matrix.col(sampler[l1])
          if i < nx-1:
            coord = matrix.col(sampler[l2])
            weight = sampler.weight(l1, coord)
            assert(abs(weight - expected) < eps)
            weight = sampler.weight(l1, ( coord + coord1 )/2.0)
            assert(abs(weight - 0.5) < eps)
          if i > 0:
            coord = matrix.col(sampler[l3])
            weight = sampler.weight(l1, coord)
            assert(abs(weight - expected) < eps)
            weight = sampler.weight(l1, ( coord1 + coord )/2.0)
            assert(abs(weight - 0.5) < eps)
          if j < ny-1:
            coord = matrix.col(sampler[l4])
            weight = sampler.weight(l1, coord)
            assert(abs(weight - expected) < eps)
            weight = sampler.weight(l1, ( coord + coord1 )/2.0)
            assert(abs(weight - 0.5) < eps)
          if j > 0:
            coord = matrix.col(sampler[l5])
            weight = sampler.weight(l1, coord)
            assert(abs(weight - expected) < eps)
            weight = sampler.weight(l1, ( coord1 + coord )/2.0)
            assert(abs(weight - 0.5) < eps)
          if k < nz-1:
            coord = matrix.col(sampler[l6])
            weight = sampler.weight(l1, coord)
            assert(abs(weight - expected) < eps)
            weight = sampler.weight(l1, ( coord + coord1 )/2.0)
            assert(abs(weight - 0.5) < eps)
          if k > 0:
            coord = matrix.col(sampler[l7])
            weight = sampler.weight(l1, coord)
            assert(abs(weight -
) < eps)
            weight = sampler.weight(l1, ( coord1 + coord )/2.0)
            assert(abs(weight - 0.5) < eps)

    print 'OK'

  def tst_self_consistent(self):
    from dials.algorithms.integration.profile import GridSampler
    width = 1000
    height = 1000
    depth = 10
    nx = 10
    ny = 10
    nz = 2
    sampler = GridSampler((width, height, depth), (nx, ny, nz))

    for i in range(len(sampler)):
      coord = sampler[i]
      index = sampler.nearest(coord)
      assert(index == i)

    print 'OK'

  def tst_pickle(self):
    from dials.algorithms.integration.profile import GridSampler
    import tempfile
    import cPickle as pickle
    width = 1000
    height = 1000
    depth = 10
    nx = 10
    ny = 10
    nz = 2
    sampler = GridSampler((width, height, depth), (nx, ny, nz))

    tf = tempfile.TemporaryFile()
    pickle.dump(sampler, tf)
    tf.flush()
    tf.seek(0)
    sampler2 = pickle.load(tf)

    assert(sampler.volume_size() == sampler2.volume_size())
    assert(sampler.grid_size() == sampler2.grid_size())

    print 'OK'

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
