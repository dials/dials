from __future__ import division
import boost.python
from dials.array_family import flex
from dials_model_serialize_ext import *

class ShoeboxExporterAux(boost.python.injector, ShoeboxFileExporter):
  ''' A class to add extra methods to the shoebox exporter class. '''

  def export(self, imageset):
    ''' Export all the shoeboxes from an imageset. '''

    from dials.util.command_line import ProgressBar

    # For each image in the sweep, get the reflections predicted to have
    # been recorded on the image and copy the pixels from the image to
    # the reflection profile image. Update a progress bar as we go along
    progress = ProgressBar(title = "Extracting shoeboxes")
    for frame, image in enumerate(imageset):

      # Ensure the image is a tuple
      if not isinstance(image, tuple):
        image = (image,)

      # Loop through all the images and add to the extractor
      for panel, im in enumerate(image):
        f, p = self.next(im)
        assert(f == frame and p == panel)

      # Update the progress bar
      progress.update(100*(frame+1) / len(imageset))

    # Finish the progress bar and return the profiles
    assert(self.finished())
    self.flush()
    progress.finished("Extracted profiles from %d frames" % len(imageset))


def extract_shoeboxes_to_file(filename, imageset, reflections):
  ''' Extract the shoeboxes to file. '''
  import cPickle as pickle

  # Get some stuff from the experiment
  detector = imageset.get_detector()
  scan = imageset.get_scan()
  frame_offset = scan.get_array_range()[0]
  num_frames = scan.get_num_images()
  assert(num_frames == len(imageset))
  num_panels = len(detector)

  # Create the shoebox file exporter
  exporter = ShoeboxFileExporter(
    filename,
    reflections['panel'],
    reflections['bbox'],
    reflections['xyzcal.px'].parts()[2],
    frame_offset,
    num_frames,
    num_panels,
    pickle.dumps(reflections, protocol=pickle.HIGHEST_PROTOCOL))

  # Write all the shoeboxes to disk
  exporter.export(imageset)
