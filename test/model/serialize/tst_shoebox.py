from __future__ import division
from dials.model.serialize.shoebox import Writer, Reader

class Test(object):

  def __init__(self):
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'FAIL: dials_regression not configured'
      exit(0)

    import os
    from dials.model.experiment.experiment_list import ExperimentListFactory
    path = os.path.join(dials_regression, 'centroid_test_data')
    self.experiments = ExperimentListFactory.from_json_file(
      os.path.join(path, "experiments.json"))

  def run(self):
    from dials.array_family import flex
    predicted = self.predict_reflections()

    filename = 'extracted.tar'
    writer = Writer(filename)
    writer.write(self.experiments[0].imageset, predicted)
    writer.close()
    del writer

    # Open the reader
    reader = Reader(filename)

    # Check the records
    maxz_list = []
    for record in reader.iter_records():
      minz, maxz = record[0]
      filename = record[1]
      maxz_list.append(maxz)
      assert(filename is not None)
    assert(len(set(maxz_list)) == len(maxz_list))
    assert(all(z1 == z2 for z1, z2 in zip(maxz_list, sorted(maxz_list))))

    # Check the shoeboxes
    all_indices = []
    for zrange, index, shoebox in reader:
      assert(len(index) == len(shoebox))
      assert(zrange[1] > zrange[0])
      assert(all(shoebox.is_consistent()))
      all_indices.extend(list(index))
    assert(len(all_indices) == len(predicted))
    assert(len(set(all_indices)) == len(all_indices))

    # Fill the expected shoeboxes
    shoeboxes = flex.basic_shoebox(predicted['panel'], predicted['bbox'])
    image_data = self.experiments[0].imageset.to_array()
    for i in range(len(shoeboxes)):
      x0, x1, y0, y1, z0, z1 = shoeboxes[i].bbox
      shoeboxes[i].data = image_data[z0:z1,y0:y1,x0:x1]

    # Check all the pixel values are as expected
    for zrange, index, shoebox in reader:
      for i, sbox1 in zip(index, shoebox):
        sbox2 = shoeboxes[i]
        assert(sbox1.data.all_eq(sbox2.data))

    # Test passed
    print 'OK'

  def predict_reflections(self):
    from dials.array_family import flex
    from math import pi
    predicted = flex.reflection_table.from_predictions(self.experiments[0])
    predicted.compute_bbox(
      self.experiments[0],
      nsigma=3,
      sigma_d = 0.058 * pi / 180.0,
      sigma_m = 0.157 * pi / 180.0)
    bbox = predicted['bbox']
    zrange = self.experiments[0].imageset.get_array_range()
    width, height = self.experiments[0].detector[0].get_image_size()
    for i in range(len(bbox)):
      x0, x1, y0, y1, z0, z1 = bbox[i]
      x0 = max(0, x0)
      x1 = min(width, x1)
      y0 = max(0, y0)
      y1 = min(height, y1)
      z0 = max(zrange[0], z0)
      z1 = min(zrange[1], z1)
      bbox[i] = (x0, x1, y0, y1, z0, z1)
    return predicted



if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
