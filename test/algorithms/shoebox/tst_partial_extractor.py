from __future__ import division

class TestSinglePanel(object):

  def __init__(self):
    ''' Create the maps and reflection list'''
    from scitbx.array_family import flex

    # Create the mask, gainmap and dark map
    xsize, ysize = (2000, 2000)
    size = xsize * ysize

    # Generate some reflections
    self.reflections = self.generate_reflections(1000)
    self.zrange = (2, 7)

  def generate_reflections(self, num):
    ''' Generate some random reflections.'''
    from dials.model.data import ReflectionList
    from random import randint
    rlist = ReflectionList(num)
    for i in range(num):
      x0 = randint(-10, 2000)
      x1 = x0 + randint(1, 10)
      y0 = randint(-10, 2000)
      y1 = y0 + randint(1, 10)
      z0 = randint(-5, 10)
      z1 = z0 + randint(1, 10)
      rlist[i].bounding_box = x0, x1, y0, y1, z0, z1
    return rlist

  def run(self):
    ''' Run the tests'''
    self.tst_initialize()
    self.tst_add_image()

  def tst_initialize(self):
    ''' Test that the indices of each pixel are the same'''
    from dials.algorithms import shoebox
    from collections import defaultdict
    from dials.array_family import flex

    # Create the extractor
    bboxes = flex.int6([r.bounding_box for r in self.reflections])
    extractor = shoebox.PartialExtractor(bboxes, self.zrange)

    assert(len(extractor.shoeboxes()) == len(extractor.shoebox_indices()))

    # Create a dictionary of flex arrays
    frames_to_reflection = defaultdict(flex.int)

    # For each reflection, Find the frames which it spans and copy an
    # index into the frame -> reflection list
    for i, r in enumerate(self.reflections):
      if r.is_valid():
        f0 = r.bounding_box[4]
        f1 = r.bounding_box[5]
        for f in range(f0, f1):
          frames_to_reflection[f].append(i)


    # Ensure both give the same indices
    for k, v in frames_to_reflection.iteritems():
      if k >= self.zrange[0] and k < self.zrange[1]:
        index1 = v
        index2 = extractor.indices(0, k)
        assert(index1 == index2)


    # Test passed
    print 'OK'

  def tst_add_image(self):
    ''' Add a load of images and then check that the values of the
    shoebox pixels are correct.'''
    from dials.algorithms import shoebox
    from dials.array_family import flex

    # Create the extractor
    bboxes = flex.int6([r.bounding_box for r in self.reflections])
    extractor = shoebox.PartialExtractor(bboxes, self.zrange)

    # Add all the images to the shoeboxes
    images = []
    for i in range(0, 10):
      image = flex.random_int_gaussian_distribution(2000*2000, 100, 10)
      image.reshape(flex.grid(2000, 2000))
      images.append(image)

      if i >= self.zrange[0] and i < self.zrange[1]:
        extractor.add_image(0, i, image)

    shoeboxes = extractor.shoeboxes()

    # Check all reflection shoeboxes have the correct values
    for r in shoeboxes:
      x0, x1, y0, y1, z0, z1 = r.bbox
      z0, z1 = r.zrange
      x00, x11, y00, y11, z00, z11 = r.bbox
      z00, z11 = r.zrange
      profile = r.data

      # Clip coordinates
      if x0 < 0: x0 = 0
      if y0 < 0: y0 = 0
      if z0 < 0: z0 = 0
      if x1 > 2000: x1 = 2000
      if y1 > 2000: y1 = 2000
      if z1 > 10: z1 = 10

      # Ensure that the values are either all zero or all the
      # same as the image at the same coordinates
      if x1 <= x0 or y1 <= y0 or z1 <= z0:
        assert(profile.all_eq(0))
      else:
        for k in range(z0, z1):
          for j in range(y0, y1):
            for i in range(x0, x1):
              kk = k - z00
              jj = j - y00
              ii = i - x00
              value1 = images[k][j, i]
              value2 = profile[kk, jj, ii]
              assert(value1 == value2)


    # Test passed
    print 'OK'

class TestMultiPanel(object):

  def __init__(self):
    ''' Create the maps and reflection list'''
    from scitbx.array_family import flex

    # Create the mask, gainmap and dark map
    xsize, ysize = (500, 500)
    size = xsize * ysize
    self.npanels = 10
    self.zrange = (2, 7)

    # Generate some reflections
    self.reflections = self.generate_reflections(1000)

  def generate_reflections(self, num):
    ''' Generate some random reflections.'''
    from dials.model.data import ReflectionList
    from random import randint
    rlist = ReflectionList(num)
    for i in range(num):
      x0 = randint(-10, 500)
      x1 = x0 + randint(1, 10)
      y0 = randint(-10, 500)
      y1 = y0 + randint(1, 10)
      z0 = randint(-5, 10)
      z1 = z0 + randint(1, 10)
      rlist[i].bounding_box = x0, x1, y0, y1, z0, z1
      rlist[i].panel = randint(0, self.npanels-1)
    return rlist

  def run(self):
    ''' Run the tests'''
    self.tst_initialize()
    self.tst_add_image()

  def tst_initialize(self):
    ''' Test that the indices of each pixel are the same'''
    from dials.algorithms import shoebox
    from collections import defaultdict
    from dials.array_family import flex

    # Create the extractor
    bboxes = flex.int6([r.bounding_box for r in self.reflections])
    panels = flex.int([r.panel_number for r in self.reflections])
    panels2 = flex.size_t(len(panels))
    for i in range(len(panels)):
      panels2[i] = panels[i]
    extractor = shoebox.PartialExtractor(panels2, bboxes,
      self.zrange, self.npanels)

    assert(len(extractor.shoeboxes()) == len(extractor.shoebox_indices()))

    # Create a dictionary of flex arrays
    frames_to_reflection = defaultdict(flex.int)

    # For each reflection, Find the frames which it spans and copy an
    # index into the frame -> reflection list
    for i, r in enumerate(self.reflections):
      if r.is_valid():
        f0 = r.bounding_box[4]
        f1 = r.bounding_box[5]
        p = r.panel_number
        for f in range(f0, f1):
          frames_to_reflection[(p, f)].append(i)

    # Ensure both give the same indices
    for k, v in frames_to_reflection.iteritems():
      if k >= self.zrange[0] and k < self.zrange[1]:
        index1 = v
        index2 = extractor.indices(k[0], k[1])
        assert(index1 == index2)


    # Test passed
    print 'OK'

  def tst_add_image(self):
    ''' Add a load of images and then check that the values of the
    shoebox pixels are correct.'''
    from dials.algorithms import shoebox
    from dials.array_family import flex

    # Create the extractor
    bboxes = flex.int6([r.bounding_box for r in self.reflections])
    panels = flex.int([r.panel_number for r in self.reflections])
    panels2 = flex.size_t(len(panels))
    for i in range(len(panels)):
      panels2[i] = panels[i]
    extractor = shoebox.PartialExtractor(panels2, bboxes,
        self.zrange, self.npanels)

    # Add all the images to the shoeboxes
    images = [[] for i in range(self.npanels)]
    for i in range(10):
      for p in range(self.npanels):
        image = flex.random_int_gaussian_distribution(500*500, 100, 10)
        image.reshape(flex.grid(500, 500))
        images[p].append(image)

        if i >= self.zrange[0] and i < self.zrange[1]:
          extractor.add_image(p, i, image)

    shoeboxes = extractor.shoeboxes()

    # Check all reflection shoeboxes have the correct values
    for r in shoeboxes:
      x0, x1, y0, y1, z0, z1 = r.bbox
      z0, z1 = r.zrange
      x00, x11, y00, y11, z00, z11 = r.bbox
      z00, z11 = r.zrange
      profile = r.data
      panel = r.panel

      # Clip coordinates
      if x0 < 0: x0 = 0
      if y0 < 0: y0 = 0
      if z0 < 0: z0 = 0
      if x1 > 500: x1 = 500
      if y1 > 500: y1 = 500
      if z1 > 10: z1 = 10

      # Ensure that the values are either all zero or all the
      # same as the image at the same coordinates
      if x1 <= x0 or y1 <= y0 or z1 <= z0:
        assert(profile.all_eq(0))
      else:
        for k in range(z0, z1):
          for j in range(y0, y1):
            for i in range(x0, x1):
              kk = k - z00
              jj = j - y00
              ii = i - x00
              value1 = images[panel][k][j, i]
              value2 = profile[kk, jj, ii]
              assert(value1 == value2)

    # Test passed
    print 'OK'

class Test:

  def __init__(self):
    self.tst_single_panel = TestSinglePanel()
    self.tst_multi_panel = TestMultiPanel()

  def run(self):
    #self.tst_single_panel.run()
    self.tst_multi_panel.run()

if __name__ == '__main__':
  test = Test()
  test.run()
