from __future__ import division
from dials.model.serialize.shoebox import Writer

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
    #blocks = self.calculate_blocks(self.sweep, 5)

    filename = 'extracted.tar'
    writer = Writer(filename)
    writer.write(self.experiments[0].imageset, predicted)
    writer.close()

    ## Open the reader
    #reader = Reader(filename)

    #all_indices = []
    #for indices, shoeboxes in reader:
      #all_indices.extend(indices)

    ## Check that each reflection is present and only once
    #all_indices = sorted(all_indices)
    #for i, j in enumerate(all_indices):
      #assert(i == j)

    # Test passed
    print 'OK'

  def calculate_blocks(self, sweep, nblocks):
    from math import ceil
    blocks = [0]
    sweep_length = len(sweep)
    block_length = int(ceil(sweep_length / nblocks))
    for i in range(nblocks):
      frame = (i + 1) * block_length
      if frame > sweep_length:
        frame = sweep_length
      blocks.append(frame)
    return blocks

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
    for i in range(len(bbox)):
      x0, x1, y0, y1, z0, z1 = bbox[i]
      z0 = max(zrange[0], z0)
      z1 = min(zrange[1], z1)
      bbox[i] = (x0, x1, y0, y1, z0, z1)
    return predicted



if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
