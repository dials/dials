from __future__ import absolute_import, division, print_function

class Test2d:

  def __init__(self):
    from dials.algorithms.image.connected_components import LabelImageStack2d
    self.size = (500, 500)
    self.label_images = LabelImageStack2d(self.size)

  def run(self):

    from scitbx.array_family import flex

    data_list = []
    mask_list = []
    for i in range(10):
      data = flex.random_int_gaussian_distribution(
          self.size[0] * self.size[1], 100, 10)
      data.reshape(flex.grid(self.size))
      mask = flex.random_bool(self.size[0] * self.size[1], 0.1)
      mask.reshape(flex.grid(self.size))
      data_list.append(data)
      mask_list.append(mask)

    for i in range(10):
      self.label_images.add_image(data_list[i], mask_list[i])

    labels = self.label_images.labels()
    coords = self.label_images.coords()
    values = self.label_images.values()

    assert(len(labels) > 0)
    assert(len(labels) == len(coords))
    assert(len(labels) == len(values))

    self.tst_coords_are_valid(mask_list, coords)
    self.tst_values_are_valid(data_list, mask_list, values)
    self.tst_labels_are_valid(data_list, mask_list, coords, labels)

  def tst_coords_are_valid(self, mask_list, coords):

    # Ensure that the values are all ok and in the right order
    vi = 0
    for k in range(10):
      ind = 0
      for j in range(self.size[0]):
        for i in range(self.size[1]):
          m = mask_list[k][ind]
          if m:
            c1 = (k, j, i)
            c2 = coords[vi]
            vi += 1
            assert(c1 == c2)
          ind += 1

    # Test passed

  def tst_values_are_valid(self, data_list, mask_list, values):

    # Ensure that the values are all ok and in the right order
    vi = 0
    for k in range(10):
      for d, m in zip(data_list[k], mask_list[k]):
        if m:
          v1 = d
          v2 = values[vi]
          vi += 1
          assert(v1 == v2)

    # Test passed


  def tst_labels_are_valid(self, data_list, mask_list, coords, labels):
    from scitbx.array_family import flex

    # Create a map of labels
    label_map = flex.int(flex.grid(10, self.size[0], self.size[1]))
    for c, l in zip(coords, labels):
      assert(c[0] >= 0 and c[0] < 10)
      assert(c[1] >= 0 and c[1] < self.size[0])
      assert(c[2] >= 0 and c[2] < self.size[1])
      label_map[c] = l

    # Ensure all labels are correct
    vi = 0
    for k in range(10):
      for j in range(self.size[0]):
        for i in range(self.size[1]):
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

class Test3d:

  def __init__(self):
    from dials.algorithms.image.connected_components import LabelImageStack3d
    self.size = (500, 500)
    self.label_images = LabelImageStack3d(self.size)

  def run(self):

    from scitbx.array_family import flex

    data_list = []
    mask_list = []
    for i in range(10):
      data = flex.random_int_gaussian_distribution(
          self.size[0] * self.size[1], 100, 10)
      data.reshape(flex.grid(self.size))
      mask = flex.random_bool(self.size[0] * self.size[1], 0.1)
      mask.reshape(flex.grid(self.size))
      data_list.append(data)
      mask_list.append(mask)

    for i in range(10):
      self.label_images.add_image(data_list[i], mask_list[i])

    labels = self.label_images.labels()
    coords = self.label_images.coords()
    values = self.label_images.values()

    assert(len(labels) > 0)
    assert(len(labels) == len(coords))
    assert(len(labels) == len(values))

    self.tst_coords_are_valid(mask_list, coords)
    self.tst_values_are_valid(data_list, mask_list, values)
    self.tst_labels_are_valid(data_list, mask_list, coords, labels)

  def tst_coords_are_valid(self, mask_list, coords):

    # Ensure that the values are all ok and in the right order
    vi = 0
    for k in range(10):
      ind = 0
      for j in range(self.size[0]):
        for i in range(self.size[1]):
          m = mask_list[k][ind]
          if m:
            c1 = (k, j, i)
            c2 = coords[vi]
            vi += 1
            assert(c1 == c2)
          ind += 1

    # Test passed

  def tst_values_are_valid(self, data_list, mask_list, values):

    # Ensure that the values are all ok and in the right order
    vi = 0
    for k in range(10):
      for d, m in zip(data_list[k], mask_list[k]):
        if m:
          v1 = d
          v2 = values[vi]
          vi += 1
          assert(v1 == v2)

    # Test passed


  def tst_labels_are_valid(self, data_list, mask_list, coords, labels):
    from scitbx.array_family import flex

    # Create a map of labels
    label_map = flex.int(flex.grid(10, self.size[0], self.size[1]))
    for c, l in zip(coords, labels):
      label_map[c] = l

    # Ensure all labels are correct
    vi = 0
    for k in range(10):
      for j in range(self.size[0]):
        for i in range(self.size[1]):
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


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test3d()
    test.run()

    test = Test2d()
    test.run()
