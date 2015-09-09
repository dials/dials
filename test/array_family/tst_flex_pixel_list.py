from __future__ import division

class Test(object):

  def __init__(self):
    pass

  def run(self):
    from dials.array_family import flex
    from dials.model.data import PixelList
    from random import randint

    size = (500, 500)
    frange1 = (10, 15)
    frange2 = (15, 20)
    n1 = 100
    n2 = 200
    v1 = flex.random_int_gaussian_distribution(n1, 100, 5).as_double()
    v2 = flex.random_int_gaussian_distribution(n2, 100, 5).as_double()
    c1 = flex.vec3_int(n1)
    c2 = flex.vec3_int(n2)

    for i in range(n1):
      c1[i] = (randint(frange1[0], frange1[1]-1),
               randint(0, size[0]-1),
               randint(0, size[1]-1))

    for i in range(n2):
      c2[i] = (randint(frange2[0], frange2[1]-1),
               randint(0, size[0]-1),
               randint(0, size[1]-1))

    fpl = flex.pixel_list(2)
    fpl[0] = PixelList(size, frange1, v1, c1)
    fpl[1] = PixelList(size, frange2, v2, c2)

    pl = fpl.merge()

    v3 = pl.values()
    c3 = pl.coords()
    assert(len(v3) == len(c3))
    assert(len(v3) == n1 + n2)
    assert(pl.first_frame() == frange1[0])
    assert(pl.last_frame() == frange2[1])

    for i in range(n1):
      assert(v3[i] == v1[i])
      assert(c3[i] == c1[i])
    for i in range(n2):
      assert(v3[n1+i] == v2[i])
      assert(c3[n1+i] == c2[i])

    print 'OK'

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
