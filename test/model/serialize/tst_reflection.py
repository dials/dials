from __future__ import division
from dials.model.serialize.reflection import Extractor

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
    predicted['ident'] = flex.size_t(range(len(predicted)))

    extractor = Extractor(self.experiments[0], predicted)

    count = 0
    for block in extractor.extract(None):
      count += len(block)
    assert(count == len(predicted))

    blocks = [(0, 3), (3, 6), (6, 9)]
    expected = []
    zcoord = predicted['xyzcal.px'].parts()[2]
    for block in blocks:
      mask = (zcoord >= block[0]) & (zcoord < block[1])
      expected.append(predicted.select(mask))

    count = 0
    for i, block in enumerate(extractor.extract(blocks)):
      z = block['xyzcal.px'].parts()[2]
      assert(all(z >= blocks[i][0]))
      assert(all(z < blocks[i][1]))
      assert(flex.size_t(sorted(block['ident'])).all_eq(expected[i]['ident']))
      count += len(block)
    assert(count == len(predicted))

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

