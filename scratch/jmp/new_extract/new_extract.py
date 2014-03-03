from __future__ import division
from contextlib import contextmanager

class ShoeboxOctreeWriter(object):

  def __init__(self, predictions, filename):
    from dials.algorithms.spatial_indexing import OctreeVec3Double

    # Create an octree of the reflection shoeboxes
    tree = OctreeVec3Double(predictions['xyzcal.px'])

    # from the octree, extract a list of the buckets
    pass

  def __del__(self):
    pass

  def add_image(self, index, image):
    pass


class ShoeboxOctreeReader(object):

  def __init__(self, filename):
    pass

  def __getitem__(self, block):
    pass


class Extractor(object):

  def __init__(self, experiment, predictions, filename=None):

    # Set the filename if there isn't one
    if filename is None:
      filename = 'extracted.tar'
    self.filename = filename

    # Extract the images to disk
    self._extract_to_disk(experiment, predictions)

    # Create the shoebox quadtree reader
    self.reader = ShoeboxOctreeReader(filename)

  def _extract_to_disk(self, experiment, predictions):
    from dials.util.command_line import ProgressBar

    # Create the shoebox quadtree writer
    writer = ShoeboxOctreeWriter(predictions, self.filename)

    # Add all the images to the writer
    progress = ProgressBar(title='Extracting shoeboxes')
    start = experiment.imageset.get_array_range()[0]
    n = len(experiment.imageset)
    for i, image in enumerate(experiment.imageset, start=start):
      writer.add_image(i, image)
      progress.update(100*(i-start+1) / n)
    progress.finished('Extracted shoeboxes')

  def __getitem__(self, block):
    return self.reader[block]

  def __call__(self, blocks):
    for block in blocks:
      yield self[block]


@contextmanager
def create_extractor(experiment, predictions, filename=None):

  # Create the extractor
  extractor = Extractor(experiment, predictions, filename)

  # Yield the extractor
  yield extractor




if __name__ == '__main__':
  from dials.model.experiment.experiment_list import ExperimentListFactory
  from dials.array_family import flex

  experiments = ExperimentListFactory.from_json_file(
    '/home/upc86896/Data/i04-BAG-training/experiments.json')

  predictions = flex.reflection_table.from_predictions(experiments)
  predictions.compute_bbox(experiments[0], nsigma=3, sigma_d=0.02, sigma_m=0.04)

  blocks = [(0, 1000, 0, 1000, 0, 540)]

  print 'Create extractor'
  with create_extractor(experiments[0], predictions) as extract:

    for block in extract(blocks):
      print block
