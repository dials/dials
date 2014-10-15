from __future__ import division
from dials.util.options import OptionParser
from dials.util.options import flatten_datablocks
from dials.util.options import flatten_experiments
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon
from scitbx.matrix import col

class Script(object):
  def __init__(self):
    # Create the parser
    self.parser = OptionParser(
      read_experiments=True,
      read_datablocks=True,
      read_datablocks_from_images=True)

  def run(self):
    params, options = self.parser.parse_args(show_diff_phil=True)
    datablocks  = flatten_datablocks(params.input.datablock)
    experiments = flatten_experiments(params.input.experiments)

    all_detectors = []
    for db in datablocks:
      all_detectors.extend(db.unique_detectors())

    all_detectors.extend(experiments.detectors())

    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')

    #for detector in all_detectors:
    for j in xrange(min(len(all_detectors),10)):
      detector = all_detectors[j]
      print detector[0].get_d_matrix()
      for i, panel in enumerate(detector):
        size = panel.get_image_size()
        p0 = col(panel.get_pixel_lab_coord((0,0)))
        p1 = col(panel.get_pixel_lab_coord((size[0]-1,0)))
        p2 = col(panel.get_pixel_lab_coord((size[0]-1,size[1]-1)))
        p3 = col(panel.get_pixel_lab_coord((0,size[1]-1)))

        v1 = p1-p0
        v2 = p3-p0
        vcen = ((v2/2) + (v1/2)) + p0

        ax.add_patch(Polygon((p0[0:2],p1[0:2],p2[0:2],p3[0:2]), closed=True, color='green', fill=False, hatch='/'))
        ax.annotate(i, vcen[0:2])

    ax.set_xlim((-100, 100))
    ax.set_ylim((-100, 100))
    plt.show()


if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
