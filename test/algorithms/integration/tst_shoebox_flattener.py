from __future__ import division


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


if __name__ == '__main__':
  test = Test()
  test.run()
