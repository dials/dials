from __future__ import absolute_import, division, print_function

class Test(object):

  def __init__(self):
    from dials.algorithms.integration.profile import GridSampler2D
    from dials.array_family import flex
    from random import randint, uniform

    # Number of reflections
    nrefl = 1000

    # Size of the images
    width = 1000
    height = 1000

    # Create the grid
    self.grid = GridSampler2D((width, height), (5,5))

    # Create the list of xyz and bboxes
    self.xyz = flex.vec3_double(nrefl)
    self.bbox = flex.int6(nrefl)
    self.panel = flex.size_t(nrefl, 0)
    for i in range(nrefl):
      x0 = randint(0,width-10)
      x1 = x0 + randint(3,10)
      y0 = randint(0,height-10)
      y1 = y0 + randint(3,10)
      z0 = randint(0,10)
      z1 = z0 + randint(1,10)
      b = x0, x1, y0, y1, z0, z1
      c = (x1 + x0) / 2, (y1 + y0) / 2, (z1 + z0) / 2
      self.xyz[i] = c
      self.bbox[i] = b

    # Create the array of shoeboxes
    self.sbox = flex.shoebox(self.panel, self.bbox)
    self.sbox.allocate()
    for i in range(len(self.sbox)):
      data = self.sbox[i].data
      for j in range(len(data)):
        data[j] = uniform(0, 100)

  def run(self):
    from dials.algorithms.integration import ShoeboxFlattener

    # Create the shoebox flattener
    flattener = ShoeboxFlattener(
      self.grid,
      self.xyz,
      self.sbox)

    # Indices in the grid
    index = flattener.index()

    # Flattened bbox
    bbox = flattener.bbox()

    # All the shoebox data
    data = [flattener.data(i) for i in range(len(flattener))]

    # Ensure that all the bboxes within a grid are the same size
    count = [0 for i in range(5*5)]
    xsize = [0 for i in range(5*5)]
    ysize = [0 for i in range(5*5)]
    for i in range(len(index)):
      j = index[i]
      b = bbox[i]
      d = data[i]
      xs = b[1] - b[0]
      ys = b[3] - b[2]
      zs = b[5] - b[4]
      ys1, xs1 = d.all()
      assert(xs == xs1)
      assert(ys == ys1)
      if count[j] == 0:
        xsize[j] = xs
        ysize[j] = ys
      else:
        assert(xsize[j] == xs)
        assert(ysize[j] == ys)

    # Check that the flattened shoebox values are correct
    for i in range(len(self.sbox)):
      b1 = self.bbox[i]
      b2 = bbox[i]
      d1 = self.sbox[i].data
      d2 = data[i]
      x0 = b2[0] - b1[0]
      x1 = b2[1] - b1[0]
      y0 = b2[2] - b1[2]
      y1 = b2[3] - b1[2]
      d11 = d1[:,y0:y1,x0:x1]
      ys1, xs1 = d2.all()
      zs2, ys2, xs2 = d11.all()
      assert(ys1 == ys2)
      assert(xs1 == xs2)
      for y in range(ys2):
        for x in range(xs2):
          value = 0
          for z in range(zs2):
            value += d11[z,y,x]
          assert(value == d2[y,x])

    # Test passed


if __name__ == '__main__':
  test = Test()
  test.run()
