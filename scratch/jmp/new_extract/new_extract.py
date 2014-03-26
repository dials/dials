from __future__ import division
from contextlib import contextmanager



#class Extractor(object):

  #def __init__(self, experiment, predictions, filename=None):

    ## Set the filename if there isn't one
    #if filename is None:
      #filename = 'extracted.tar'
    #self.filename = filename

    ## Extract the images to disk
    #self._extract_to_disk(experiment, predictions)

    ## Create the shoebox quadtree reader
    #self.reader = ShoeboxOctreeReader(filename)

  #def _extract_to_disk(self, experiment, predictions):
    #from dials.util.command_line import ProgressBar

    ## Create the shoebox quadtree writer
    #writer = ShoeboxOctreeWriter(predictions, self.filename)

    ## Add all the images to the writer
    #progress = ProgressBar(title='Extracting shoeboxes')
    #start = experiment.imageset.get_array_range()[0]
    #n = len(experiment.imageset)
    #for i, image in enumerate(experiment.imageset, start=start):
      #writer.add_image(i, image)
      #progress.update(100*(i-start+1) / n)
    #progress.finished('Extracted shoeboxes')

  #def __getitem__(self, block):
    #return self.reader[block]

  #def __call__(self, blocks):
    #for block in blocks:
      #yield self[block]


#@contextmanager
#def create_extractor(experiment, predictions, filename=None):

  ## Create the extractor
  #extractor = Extractor(experiment, predictions, filename)

  ## Yield the extractor
  #yield extractor


class ShoeboxWriter(object):

  def __init__(self, filename):
    self.filename = filename

  def write(self, predictions, imageset):
    from dials.util.command_line import ProgressBar
    zrange = imageset.get_array_range()
    p = ProgressBar(title="Extracting shoeboxes")
    for z, image in enumerate(imageset, start=zrange[0]):
      batch = self._add_image(z, image)
      self._write_pickle(batch)
      p.update(100.0*(z - zrange[0])/len(imageset))
    p.finish("Extracted shoeboxes")
    self._write_predictions(predictions)

  def _add_image(self, z, image):
    pass

  def _write_pickle(self, batch):
    pass

  def _write_predictions(self, predictions):
    pass


@contextmanager
def open_shoebox_writer(filename):
  writer = ShoeboxWriter(filename)
  yield writer


if __name__ == '__main__':
  from dials.model.experiment.experiment_list import ExperimentListFactory
  from dials.array_family import flex

  experiments = ExperimentListFactory.from_json_file(
    '/home/upc86896/Data/Data/i04-BAG-training/dials_processed/experiments.json')

  predictions = flex.reflection_table.from_predictions(experiments[0])
  predictions.compute_bbox(experiments[0], nsigma=3, sigma_d=0.024,
                           sigma_m=0.044)

  zeta = predictions.compute_zeta(experiments[0])
  mask = flex.abs(zeta) < 0.05
  predictions.del_selected(mask)

  with open_shoebox_writer("extracted.tar") as writer:
    writer.write(predictions, experiments[0].imageset)
