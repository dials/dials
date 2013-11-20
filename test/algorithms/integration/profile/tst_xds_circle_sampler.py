from __future__ import division

class Test(object):

  def __init__(self):
    pass

  def run(self):
    self.tst_getters()
    self.tst_detector_area()
    self.tst_indexing()
    self.tst_nearest()
    self.tst_nearest_n()
    self.tst_self_consistent()
    self.tst_z_index()
    self.tst_pickle()

  def tst_getters(self):
    from math import sqrt
    from dials.algorithms.integration.profile import XdsCircleSampler
    width = 1000
    height = 1000
    depth = 10
    nz = 2
    sampler = XdsCircleSampler((width, height, depth), nz)
    volume_size = sampler.volume_size()
    image_centre = sampler.image_centre()
    r0 = sampler.r0()
    r1 = sampler.r1()
    r2 = sampler.r2()
    size = len(sampler)

    assert(width == volume_size[0])
    assert(height == volume_size[1])
    assert(depth == volume_size[2])
    assert(width // 2 == image_centre[0])
    assert(height // 2 == image_centre[1])
    assert(r0 == min([width // 2, height // 2]))
    assert(r1 == r0 / 3.0)
    assert(r2 == r1 * sqrt(5.0))
    assert(9 * nz == size)
    print 'OK'

  def tst_detector_area(self):
    from dials.algorithms.integration.profile import XdsCircleSampler
    from math import sin, cos, pi
    from scitbx.array_family import flex
    width = 1000
    height = 1000
    depth = 10
    nz = 2
    sampler = XdsCircleSampler((width, height, depth), nz)
    im = flex.int(flex.grid(height, width))
    for j in range(height):
      for i in range(width):
        im[j,i] = sampler.nearest((i, j, 0))

    assert(im[height//2, width//2] == 0)
    assert(im[height//2, width-1] == 1)
    assert(im[height-1, width-1] == 2)
    assert(im[height-1, width//2] == 3)
    assert(im[height-1, 0] == 4)
    assert(im[height//2, 0] == 5)
    assert(im[0, 0] == 6)
    assert(im[0, width//2] == 7)
    assert(im[0, width-1] == 8)
    print 'OK'
#        from matplotlib import pylab
#        pylab.imshow(im.as_numpy_array())
#        pylab.show()


  def tst_indexing(self):
    from dials.algorithms.integration.profile import XdsCircleSampler
    from math import sin, cos, pi
    width = 1000
    height = 1000
    depth = 10
    nz = 2
    sampler = XdsCircleSampler((width, height, depth), nz)
    zstep = depth / nz
    eps = 1e-10

    dt = 2 * pi / 8
    r2 = sampler.r2()

    xp = [width / 2.0] + [width / 2.0 + r2 * cos(i * dt) for i in range(8)]
    yp = [height / 2.0] + [height / 2.0 + r2 * sin(i * dt) for i in range(8)]
    zind = [[k] * 9 for k in range(nz)]
    zind = [i for j in zind for i in j]
    zp = [(z + 0.5) * zstep for z in zind]

    for x0, y0, z0, (x1, y1, z1) in zip(xp, yp, zp, sampler):
      assert(abs(x0 - x1) <= eps)
      assert(abs(y0 - y1) <= eps)
      assert(abs(z0 - z1) <= eps)

    print 'OK'


  def tst_nearest(self):
    from math import sqrt, atan2, pi, floor
    from random import randint
    from dials.algorithms.integration.profile import XdsCircleSampler
    width = 1000
    height = 1000
    depth = 10
    nz = 2
    sampler = XdsCircleSampler((width, height, depth), nz)
    xc, yc = sampler.image_centre()
    r1 = sampler.r1()

    for i in range(1000):
      x = randint(0, 1000)
      y = randint(0, 1000)
      z = randint(0, 10)

      r = sqrt((x - xc)**2 + (y - yc)**2)
      if r < r1:
        index00 = 0
      else:
        t = atan2(y - yc, x - xc)
        ai = int(floor(t * 8 / (2 * pi) + 0.5)) % 8
        index00 = ai + 1

      index01 = int(z * nz / 10)
      if index01 < 0:
        index01 = 0
      if index01 >= 2:
        index01 = 1

      index0 = index00 + index01 * 9
      index1 = sampler.nearest((x, y, z))
      assert(index0 == index1)

    print 'OK'

  def tst_nearest_n(self):
    from math import sqrt, atan2, pi, floor
    from random import randint
    from dials.algorithms.integration.profile import XdsCircleSampler
    width = 1000
    height = 1000
    depth = 10
    nz = 2
    sampler = XdsCircleSampler((width, height, depth), nz)
    xc, yc = sampler.image_centre()
    r1 = sampler.r1()

    for i in range(1000):
      x = randint(0, 1000)
      y = randint(0, 1000)
      z = randint(0, 10)

      r = sqrt((x - xc)**2 + (y - yc)**2)
      if r < r1:
        index00 = 0
      else:
        t = atan2(y - yc, x - xc)
        ai = int(floor(t * 8 / (2 * pi) + 0.5)) % 8
        index00 = ai + 1

      index01 = int(z * nz / 10)
      if index01 < 0:
        index01 = 0
      if index01 >= 2:
        index01 = 1

      index0 = index00 + index01 * 9
      index1 = sampler.nearest_n((x, y, z))
      assert(index0 == index1[0])
      if index0 % 9 == 0:
        assert(len(index1) == 9)
        assert(all([idx == index0 + i for i, idx in enumerate(index1)]))
      else:
        assert(len(index1) == 4)
        assert(index1[1] == (index0 // 9) * 9)
        if (index0 % 9) == 1:
          assert(index1[2] == index0 + 1)
          assert(index1[3] == index0 + 7)
        elif (index0 % 9) == 8:
          assert(index1[2] == index0 - 7)
          assert(index1[3] == index0 - 1)
        else:
          assert(index1[2] == index0 + 1)
          assert(index1[3] == index0 - 1)

    print 'OK'

  def tst_self_consistent(self):
    from dials.algorithms.integration.profile import XdsCircleSampler
    width = 1000
    height = 1000
    depth = 10
    nz = 2
    sampler = XdsCircleSampler((width, height, depth), nz)

    for i in range(len(sampler)):
      coord = sampler[i]
      index = sampler.nearest(coord)
      assert(index == i)

    print 'OK'

  def tst_z_index(self):
    from math import sqrt, atan2, pi, floor
    from random import randint
    from dials.algorithms.integration.profile import XdsCircleSampler
    width = 1000
    height = 1000
    depth = 10
    nz = 2
    sampler = XdsCircleSampler((width, height, depth), nz)
    assert((sampler.nearest((500, 500, 0.0)) / 9) == 0)
    assert((sampler.nearest((500, 500, 1.0)) / 9) == 0)
    assert((sampler.nearest((500, 500, 2.0)) / 9) == 0)
    assert((sampler.nearest((500, 500, 3.0)) / 9) == 0)
    assert((sampler.nearest((500, 500, 4.0)) / 9) == 0)
    assert((sampler.nearest((500, 500, 4.5)) / 9) == 0)
    assert((sampler.nearest((500, 500, 5.0)) / 9) == 1)
    assert((sampler.nearest((500, 500, 5.5)) / 9) == 1)
    assert((sampler.nearest((500, 500, 6.0)) / 9) == 1)
    assert((sampler.nearest((500, 500, 7.0)) / 9) == 1)
    assert((sampler.nearest((500, 500, 8.0)) / 9) == 1)
    assert((sampler.nearest((500, 500, 9.0)) / 9) == 1)
    assert((sampler.nearest((500, 500, 10.0)) / 9) == 1)

    print 'OK'

  def tst_pickle(self):
    from dials.algorithms.integration.profile import XdsCircleSampler
    import tempfile
    import cPickle as pickle
    width = 1000
    height = 1000
    depth = 10
    nz = 2
    sampler = XdsCircleSampler((width, height, depth), nz)

    tf = tempfile.TemporaryFile()
    pickle.dump(sampler, tf)
    tf.flush()
    tf.seek(0)
    sampler2 = pickle.load(tf)

    assert(sampler.volume_size() == sampler2.volume_size())
    assert(sampler.num_z() == sampler2.num_z())

    print 'OK'

if __name__ == '__main__':
  test = Test()
  test.run()
