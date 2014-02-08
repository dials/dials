from __future__ import division
from dials.model.serialize.partial_shoebox import Reader, Writer

class Test(object):

  def __init__(self):
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'FAIL: dials_regression not configured'
      exit(0)

    import os
    from dials.model.serialize import load, dump
    from cctbx.crystal.crystal_model.serialize import load_crystal
    path = os.path.join(dials_regression, 'centroid_test_data')
    self.sweep = load.sweep(os.path.join(path, 'sweep.json'))
    self.crystal = load_crystal(os.path.join(path, 'crystal.json'))

  def run(self):
    from dials.algorithms.integration import PartialProfileExtractor
    from dials.array_family import flex

    predicted = self.predict_reflections()
    blocks = self.calculate_blocks(self.sweep, 5)

    filename = 'extracted.tar'
    writer = Writer(filename, predicted, blocks)

    panels = flex.size_t([p.panel_number for p in predicted])
    bboxes = flex.int6([p.bounding_box for p in predicted])

    # Extract the frames
    for i, (b0, b1) in enumerate(zip(blocks[:-1], blocks[1:])):
      extract = PartialProfileExtractor(self.sweep[b0:b1])
      writer.write(i, *extract(panels, bboxes))
    writer.close()

    # Open the reader
    reader = Reader(filename)

    all_indices = []
    for indices, shoeboxes in reader:
      all_indices.extend(indices)

    # Check that each reflection is present and only once
    all_indices = sorted(all_indices)
    for i, j in enumerate(all_indices):
      assert(i == j)

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
    from dials.util.command_line import Command
    from dials.algorithms.integration import ReflectionPredictor
    from dials.algorithms.shoebox import BBoxCalculator
    from dials.model.data import ReflectionList
    from dials.array_family import flex

    predict = ReflectionPredictor()
    predicted = predict(self.sweep, self.crystal)

    n_sigma = 3

    # Create the bbox calculator
    compute_bbox = BBoxCalculator(
        self.sweep.get_beam(), self.sweep.get_detector(),
        self.sweep.get_goniometer(), self.sweep.get_scan(),
        n_sigma * self.sweep.get_beam().get_sigma_divergence(deg=False),
        n_sigma * self.crystal.get_mosaicity(deg=False))

    # Calculate the bounding boxes of all the reflections
    Command.start('Calculating bounding boxes')
    compute_bbox(predicted)
    Command.end('Calculated {0} bounding boxes'.format(len(predicted)))

    return predicted



if __name__ == '__main__':
  test = Test()
  test.run()
