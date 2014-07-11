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
    from dials.algorithms.integration.profile import GridSampler2D
    width = 1000
    height = 1000
    nx = 5
    ny = 5
    sampler = GridSampler2D((width, height), (nx, ny))
    image_size = sampler.image_size()
    grid_size = sampler.grid_size()
    step_size = sampler.step_size()
    size = len(sampler)

    assert(width == image_size[0])
    assert(height == image_size[1])
    assert(nx == grid_size[0])
    assert(ny == grid_size[1])
    assert(step_size[0] == width / nx)
    assert(step_size[1] == height / ny)
    assert(nx * ny == size)
    print 'OK'

  def tst_indexing(self):
    from dials.algorithms.integration.profile import GridSampler2D
    width = 1000
    height = 1000
    nx = 5
    ny = 5
    sampler = GridSampler2D((width, height), (nx, ny))
    xstep, ystep = sampler.step_size()
    xind = [[i for i in range(nx)]] * ny
    yind = [[j] * nx for j in range(ny)]
    xind = [i for j in xind for i in j]
    yind = [i for j in yind for i in j]

    xp = [(x + 0.5) * xstep for x in xind]
    yp = [(y + 0.5) * ystep for y in yind]

    eps = 1e-10

    for x0, y0, (x1, y1) in zip(xp, yp, sampler):
      assert(abs(x0 - x1) <= eps)
      assert(abs(y0 - y1) <= eps)

    print 'OK'


  def tst_nearest(self):
    from random import randint
    from dials.algorithms.integration.profile import GridSampler2D
    width = 1000
    height = 1000
    nx = 5
    ny = 5
    sampler = GridSampler2D((width, height), (nx, ny))

    for i in range(1000):
      x = randint(0, 1000-1)
      y = randint(0, 1000-1)
      z = randint(0, 10)
      i = int((x+0.5) * nx // 1000)
      j = int((y+0.5) * ny // 1000)
      if i >= nx:
        i = nx - 1
      if j >= ny:
        j = ny - 1
      index0 = i + j * nx
      index1 = sampler.nearest((x, y))
      assert(index0 == index1)

    print 'OK'

  def tst_nearest_n(self):
    from random import randint, uniform
    from dials.algorithms.integration.profile import GridSampler2D
    width = 1000
    height = 1000
    nx = 5
    ny = 5
    sampler = GridSampler2D((width, height), (nx, ny))

    def first_quad():
      i = randint(0, nx-1)
      j = randint(0, ny-1)
      x = (i + uniform(0, 0.5)) * 1000.0 / nx
      y = (j + uniform(0, 0.5)) * 1000.0 / ny
      z = randint(0, 10)
      indices = [i+j*nx]
      if i > 0:
        indices.append((i-1)+j*nx)
      if j > 0:
        indices.append(i+(j-1)*nx)
      if j > 0 and i > 0:
        indices.append((i-1)+(j-1)*nx)
      return x, y, z, indices

    def second_quad():
      i = randint(0, nx-1)
      j = randint(0, ny-1)
      x = (i + uniform(0, 0.5) + 0.5) * 1000.0 / nx
      y = (j + uniform(0, 0.5)) * 1000.0 / ny
      z = randint(0, 10)
      indices = [i+j*nx]
      if i < nx-1:
        indices.append((i+1)+j*nx)
      if j > 0:
        indices.append(i+(j-1)*nx)
      if j > 0 and i < nx-1:
        indices.append((i+1)+(j-1)*nx)
      return x, y, z, indices

    def third_quad():
      i = randint(0, nx-1)
      j = randint(0, ny-1)
      x = (i + uniform(0, 0.5)) * 1000.0 / nx
      y = (j + uniform(0, 0.5) + 0.5) * 1000.0 / ny
      z = randint(0, 10)
      indices = [i+j*nx]
      if i > 0:
        indices.append((i-1)+j*nx)
      if j < ny -1:
        indices.append(i+(j+1)*nx)
      if i > 0 and j < ny-1:
        indices.append((i-1)+(j+1)*nx)
      return x, y, z, indices

    def fourth_quad():
      i = randint(0, nx-1)
      j = randint(0, ny-1)
      x = (i + uniform(0, 0.5) + 0.5) * 1000.0 / nx
      y = (j + uniform(0, 0.5) + 0.5) * 1000.0 / ny
      z = randint(0, 10)
      indices = [i+j*nx]
      if i < nx - 1:
        indices.append((i+1)+j*nx)
      if j < ny -1:
        indices.append(i+(j+1)*nx)
      if i < nx -1 and j < ny-1:
        indices.append((i+1)+(j+1)*nx)
      return x, y, z, indices

    for i in range(250):
      x, y, z, indices0 = first_quad()
      indices1 = sampler.nearest_n((x, y))
      assert(len(indices1) == len(indices0))
      for j1, j2 in zip(sorted(indices1), sorted(indices0)):
        assert(j1 == j2)
      x, y, z, indices0 = second_quad()
      indices1 = sampler.nearest_n((x, y))
      assert(len(indices1) == len(indices0))
      for j1, j2 in zip(sorted(indices1), sorted(indices0)):
        assert(j1 == j2)
      x, y, z, indices0 = third_quad()
      indices1 = sampler.nearest_n((x, y))
      assert(len(indices1) == len(indices0))
      for j1, j2 in zip(sorted(indices1), sorted(indices0)):
        assert(j1 == j2)
      x, y, z, indices0 = fourth_quad()
      indices1 = sampler.nearest_n((x, y))
      assert(len(indices1) == len(indices0))
      for j1, j2 in zip(sorted(indices1), sorted(indices0)):
        assert(j1 == j2)

    print 'OK'

  def tst_weights(self):
    from dials.algorithms.integration.profile import GridSampler2D
    from scitbx import matrix
    from math import log, exp
    width = 1000
    height = 1000
    nx = 5
    ny = 5
    sampler = GridSampler2D((width, height), (nx, ny))

    # Check the weight at the coord in 1.0
    eps = 1e-7
    for i in range(len(sampler)):
      coord = sampler[i]
      weight = sampler.weight(i, coord)
      assert(abs(weight - 1.0) < eps)

    # Ensure we get the expected weight at the next grid point at half way
    # between grid points
    expected = exp(-4.0 * log(2.0))
    for j in range(ny):
      for i in range(nx):
        l1 = (i + 0) + (j + 0) * nx
        l2 = (i + 1) + (j + 0) * nx
        l3 = (i - 1) + (j + 0) * nx
        l4 = (i + 0) + (j + 1) * nx
        l5 = (i + 0) + (j - 1) * nx
        l6 = (i + 0) + (j + 0) * nx
        l7 = (i + 0) + (j + 0) * nx
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

    print 'OK'

  def tst_self_consistent(self):
    from dials.algorithms.integration.profile import GridSampler2D
    width = 1000
    height = 1000
    nx = 5
    ny = 5
    sampler = GridSampler2D((width, height), (nx, ny))

    for i in range(len(sampler)):
      coord = sampler[i]
      index = sampler.nearest(coord)
      assert(index == i)

    print 'OK'

  def tst_pickle(self):
    from dials.algorithms.integration.profile import GridSampler2D
    import tempfile
    import cPickle as pickle
    width = 1000
    height = 1000
    nx = 5
    ny = 5
    sampler = GridSampler2D((width, height), (nx, ny))

    tf = tempfile.TemporaryFile()
    pickle.dump(sampler, tf)
    tf.flush()
    tf.seek(0)
    sampler2 = pickle.load(tf)

    assert(sampler.image_size() == sampler2.image_size())
    assert(sampler.grid_size() == sampler2.grid_size())

    print 'OK'

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
