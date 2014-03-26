
from __future__ import division

class Test(object):

  def __init__(self):
    from scitbx.array_family import flex

    self.width = 500
    self.height = 600
    self.zrange = (-5, 5)

    self.images = []
    wh = self.width * self.height
    for z in range(*self.zrange):
      image = flex.int(flex.grid(self.height, self.width))
      for j in range(self.height):
        for i in range(self.width):
          image[j,i] = i + j * self.width + (z - self.zrange[0]) * wh
      self.images.append(image)


  def run(self):
    from dials.model.serialize import BlockList
    from dials.array_family import flex
    from random import randint
    from collections import defaultdict
    panel = flex.size_t(100, 0)
    bbox = flex.int6(100)
    for i in range(100):
      x0 = randint(0, self.width-10)
      x1 = x0 + randint(1, 10)
      y0 = randint(0, self.height-10)
      y1 = y0 + randint(1, 10)
      z0 = randint(self.zrange[0], self.zrange[1] - 5)
      z1 = z0 + randint(1, 5)
      bbox[i] = (x0, x1, y0, y1, z0, z1)

    expected = defaultdict(list)
    for i in range(100):
      expected[bbox[i][5]].append(i)

    blist = BlockList(panel, bbox, self.zrange)
    for z, image in enumerate(self.images, start=self.zrange[0]):
      zrange, indices, shoeboxes = blist.next(image)
      exp = expected[z+1]
      assert(len(indices) == len(shoeboxes))
      assert(len(exp) == len(indices))
      assert(all(i1 == i2 for i1, i2 in zip(exp, sorted(indices))))
      minz, maxz = self.zrange[0], z+1
      for i in expected[z+1]:
        if bbox[i][4] < minz:
          minz = bbox[i][4]
      assert(minz == zrange[0] and maxz == zrange[1])

      wh = self.width * self.height
      for sbox in shoeboxes:
        data = sbox.data
        x0, x1, y0, y1, z0, z1 = sbox.bbox
        assert(sbox.is_consistent())
        for k in range(data.all()[0]):
          for j in range(data.all()[1]):
            for i in range(data.all()[2]):
              value = (i+x0)+(j+y0)*self.width+(k+z0-self.zrange[0])*wh
              assert(value == data[k,j,i])
    assert(blist.empty())

    print 'OK'

if __name__ == '__main__':
  test = Test()
  test.run()
