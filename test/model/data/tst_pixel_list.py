from __future__ import division

class Test(object):

  def __init__(self):
    pass

  def run(self):
    self.tst_construct()
    self.tst_add_image()
    self.tst_labels_3d()
    self.tst_labels_2d()
    self.tst_with_no_points()

  def tst_construct(self):

    from dials.model.data import PixelList
    from scitbx.array_family import flex
    from random import randint
    size = (2000, 2000)
    sf = 10
    pl = PixelList(size, sf)
    assert(pl.size() == size)
    assert(pl.first_frame() == sf)
    assert(pl.last_frame() == sf)
    assert(pl.num_frames() == 0)
    assert(pl.frame_range() == (sf, sf))

    frame_range = (10, 20)
    values = flex.double(range(100))
    coords = flex.vec3_int(100)
    for i in range(100):
      coords[i] = (
          randint(10, 20-1),
          randint(0, 2000-1),
          randint(0, 2000-1))

    pl = PixelList(size, frame_range, values, coords)
    assert(pl.size() == size)
    assert(pl.first_frame() == frame_range[0])
    assert(pl.last_frame() == frame_range[1])
    assert(pl.num_frames() == frame_range[1] - frame_range[0])
    assert(pl.frame_range() == frame_range)
    assert(len(pl.values()) == 100)
    assert(len(pl.coords()) == 100)
    for i in range(100):
      assert(pl.values()[i] == values[i])
      assert(pl.coords()[i] == coords[i])
    print 'OK'

  def tst_add_image(self):
    from dials.model.data import PixelList
    from scitbx.array_family import flex
    size = (2000, 2000)
    sf = 10
    pl = PixelList(size, sf)

    count = 0
    for i in range(3):
      image = flex.random_int_gaussian_distribution(size[0]*size[1], 100, 5)
      mask = flex.random_bool(size[0]*size[1], 0.5)
      image.reshape(flex.grid(size))
      mask.reshape(flex.grid(size))
      count += len(mask.as_1d().select(mask.as_1d()))
      pl.add_image(image, mask)
    assert(len(pl.values()) == count)

    print 'OK'

  def tst_labels_3d(self):
    from dials.model.data import PixelList
    from scitbx.array_family import flex
    size = (500, 500)
    sf = 0
    pl = PixelList(size, sf)

    count = 0
    mask_list = []
    for i in range(3):
      image = flex.random_int_gaussian_distribution(size[0]*size[1], 100, 5)
      mask = flex.random_bool(size[0]*size[1], 0.5)
      image.reshape(flex.grid(size))
      mask.reshape(flex.grid(size))
      count += len(mask.as_1d().select(mask.as_1d()))
      pl.add_image(image, mask)
      mask_list.append(mask)

    coords = pl.coords()
    labels = pl.labels_3d()

    # Create a map of labels
    label_map = flex.int(flex.grid(3, size[0], size[1]))
    for c, l in zip(coords, labels):
      label_map[c] = l

    # Ensure all labels are correct
    vi = 0
    for k in range(3):
      for j in range(size[0]):
        for i in range(size[1]):
          if mask_list[k][j,i]:

            l1 = labels[vi]
            if k > 0 and mask_list[k-1][j,i]:
              l2 = label_map[k-1,j,i]
              assert(l2 == l1)
            if j > 0 and mask_list[k][j-1,i]:
              l2 = label_map[k,j-1,i]
              assert(l2 == l1)
            if i > 0 and mask_list[k][j,i-1]:
              l2 = label_map[k,j,i-1]
              assert(l2 == l1)
            vi += 1

    # Test passed
    print 'OK'

  def tst_labels_2d(self):
    from dials.model.data import PixelList
    from scitbx.array_family import flex
    size = (500, 500)
    sf = 0
    pl = PixelList(size, sf)

    count = 0
    mask_list = []
    for i in range(3):
      image = flex.random_int_gaussian_distribution(size[0]*size[1], 100, 5)
      mask = flex.random_bool(size[0]*size[1], 0.5)
      image.reshape(flex.grid(size))
      mask.reshape(flex.grid(size))
      count += len(mask.as_1d().select(mask.as_1d()))
      pl.add_image(image, mask)
      mask_list.append(mask)

    coords = pl.coords()
    labels = pl.labels_2d()

    # Create a map of labels
    label_map = flex.int(flex.grid(3, size[0], size[1]))
    for c, l in zip(coords, labels):
      label_map[c] = l

    # Ensure all labels are correct
    vi = 0
    for k in range(3):
      for j in range(size[0]):
        for i in range(size[1]):
          if mask_list[k][j,i]:

            l1 = labels[vi]
            if k > 0 and mask_list[k-1][j,i]:
              l2 = label_map[k-1,j,i]
              assert(l2 != l1)
            if j > 0 and mask_list[k][j-1,i]:
              l2 = label_map[k,j-1,i]
              assert(l2 == l1)
            if i > 0 and mask_list[k][j,i-1]:
              l2 = label_map[k,j,i-1]
              assert(l2 == l1)
            vi += 1

    # Test passed
    print 'OK'

  def tst_with_no_points(self):

    from dials.model.data import PixelList
    from scitbx.array_family import flex
    size = (500, 500)
    sf = 0
    pl = PixelList(size, sf)

    count = 0
    mask_list = []
    for i in range(3):
      image = flex.random_int_gaussian_distribution(size[0]*size[1], 100, 5)
      mask = flex.bool(size[0]*size[0], False)
      image.reshape(flex.grid(size))
      mask.reshape(flex.grid(size))
      count += len(mask.as_1d().select(mask.as_1d()))
      pl.add_image(image, mask)
      mask_list.append(mask)

    coords = pl.coords()
    labels1 = pl.labels_2d()
    labels2 = pl.labels_2d()

    assert len(coords) == 0
    assert len(labels1) == 0
    assert len(labels2) == 0

    print 'OK'

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
