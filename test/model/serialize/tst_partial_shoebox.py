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
    from dials.algorithms.shoebox import PartialProfileExtractor
    from dials.array_family import flex
    from dials.framework.registry import Registry
    registry = Registry()
    params = registry.params()
    params.integration.shoebox.sigma_b = 0.058
    params.integration.shoebox.sigma_m = 0.157

    predicted = self.predict_reflections()
    blocks = self.calculate_blocks(self.sweep, 5)

    filename = 'extracted.tar'
    writer = Writer(filename, predicted, blocks)

    panels = predicted['panel']
    bboxes = predicted['bbox']

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
    from dials.array_family import flex
    from dials.model.experiment.experiment_list import Experiment
    from dials.model.experiment.experiment_list import ExperimentList

    exlist = ExperimentList()
    exlist.append(Experiment(
      imageset=self.sweep,
      beam=self.sweep.get_beam(),
      detector=self.sweep.get_detector(),
      goniometer=self.sweep.get_goniometer(),
      scan=self.sweep.get_scan(),
      crystal=self.crystal))

    predicted = flex.reflection_table.from_predictions(exlist[0])
    predicted.compute_bbox(exlist[0], nsigma=3)
    return predicted



if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
