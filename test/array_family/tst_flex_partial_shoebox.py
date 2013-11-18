from __future__ import division

class Test(object):

  def __init__(self):
    pass

  def run(self):
    self.tst_consistent()
    self.tst_bounding_boxes()
    self.tst_zranges()
    self.tst_merge()

  def tst_consistent(self):

    from random import randint
    from dials.model.data import PartialShoebox
    from dials.array_family import flex

    shoebox = flex.partial_shoebox(10)

    for i in range(10):
      x0 = randint(0, 1000)
      y0 = randint(0, 1000)
      z0 = randint(0, 1000)
      x1 = randint(1, 10) + x0
      y1 = randint(1, 10) + y0
      z1 = randint(1, 10) + z0

      shoebox[i] = PartialShoebox((x0, x1, y0, y1, z0, z1), (z0, z1))

    assert(shoebox.is_consistent() == flex.bool(10, False))

    for i in range(10):
      shoebox[i].allocate()

    assert(shoebox.is_consistent() == flex.bool(10, True))

    for i in [0, 2, 4, 6, 8]:
      shoebox[i].data.resize(flex.grid(10, 10, 10))

    assert(shoebox.is_consistent() == flex.bool([False, True] * 5))

    for i in range(10):
      shoebox[i].deallocate()

    assert(shoebox.is_consistent() == flex.bool(10, False))

    # Test passed
    print 'OK'

  def tst_bounding_boxes(self):
    from dials.model.data import PartialShoebox
    from random import randint, sample
    from dials.array_family import flex

    shoebox = flex.partial_shoebox(10)
    bbox = flex.int6(10)
    for i in range(10):
      x0 = randint(0, 90)
      y0 = randint(0, 90)
      z0 = randint(0, 90)
      x1 = randint(1, 10) + x0
      y1 = randint(1, 10) + y0
      z1 = randint(1, 10) + z0
      bbox[i] = (x0, x1, y0, y1, z0, z1)
      shoebox[i] = PartialShoebox(bbox[i], (z0, z1))

    bbox2 = shoebox.bounding_boxes()
    for i in range(10):
      assert(bbox2[i] == bbox[i])

    # Test passed
    print 'OK'

  def tst_zranges(self):
    from dials.model.data import PartialShoebox
    from random import randint, sample
    from dials.array_family import flex, shared

    shoebox = flex.partial_shoebox(10)
    bbox = flex.int6(10)
    zrange = shared.tiny_int_2(10)
    for i in range(10):
      x0 = randint(0, 90)
      y0 = randint(0, 90)
      z0 = randint(0, 90)
      x1 = randint(1, 10) + x0
      y1 = randint(1, 10) + y0
      z1 = randint(1, 10) + z0
      bbox[i] = (x0, x1, y0, y1, z0, z1)
      zrange[i] = (z0, z1)
      shoebox[i] = PartialShoebox(bbox[i], (z0, z1))

    zrange2 = shoebox.zranges()
    for i in range(10):
      assert(zrange2[i] == zrange[i])

    # Test passed
    print 'OK'

  def tst_merge(self):
    from dials.model.data import PartialShoebox
    from random import randint, sample
    from dials.array_family import flex

    # This should work
    shoeboxes = flex.partial_shoebox(3)
    bbox = (10, 20, 10, 20, 0, 9)
    shoeboxes[0] = PartialShoebox(bbox, (0, 3))
    shoeboxes[1] = PartialShoebox(bbox, (3, 6))
    shoeboxes[2] = PartialShoebox(bbox, (6, 9))
    shoeboxes[0].allocate()
    shoeboxes[1].allocate()
    shoeboxes[2].allocate()
    shoeboxes[0].data = flex.int(flex.grid(3, 10, 10), 1)
    shoeboxes[1].data = flex.int(flex.grid(3, 10, 10), 2)
    shoeboxes[2].data = flex.int(flex.grid(3, 10, 10), 3)
    shoebox = shoeboxes.merge((0, 9))
    assert(shoebox.is_consistent())
    assert(shoebox.size() == (9, 10, 10))
    assert(shoebox.data[0:3,:,:].all_eq(1))
    assert(shoebox.data[3:6,:,:].all_eq(2))
    assert(shoebox.data[6:9,:,:].all_eq(3))

    # This shouldn't work
    shoeboxes = flex.partial_shoebox(3)
    bbox = (10, 20, 10, 20, 0, 9)
    shoeboxes[0] = PartialShoebox(bbox, (0, 3))
    shoeboxes[1] = PartialShoebox(bbox, (2, 6))
    shoeboxes[2] = PartialShoebox(bbox, (6, 9))
    shoeboxes[0].allocate()
    shoeboxes[1].allocate()
    shoeboxes[2].allocate()
    shoeboxes[0].data = flex.int(flex.grid(3, 10, 10), 1)
    shoeboxes[1].data = flex.int(flex.grid(3, 10, 10), 2)
    shoeboxes[2].data = flex.int(flex.grid(3, 10, 10), 3)

    try:
      shoebox = shoeboxes.merge((0, 9))
      assert(False)
    except Exception:
      pass

    print 'OK'

if __name__ == '__main__':
  test = Test()
  test.run()
