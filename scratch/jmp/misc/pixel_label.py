


if __name__ == '__main__':

  from dials.model.experiment.experiment_list import ExperimentListFactory
  path = "/home/upc86896/Projects/cctbx/sources/dials_regression/centroid_test_data/experiments.json"
  exlist = ExperimentListFactory.from_json_file(path)

  # The models
  detector = exlist[0].detector
  beam = exlist[0].beam
  gonio = exlist[0].goniometer
  scan = exlist[0].scan
  crystal = exlist[0].crystal

  from dials.array_family import flex
  from scitbx import matrix
  from math import floor
  width, height = detector[0].get_image_size()
  H = flex.int(flex.grid(height,width))
  K = flex.int(flex.grid(height,width))
  L = flex.int(flex.grid(height,width))

  A = crystal.get_A()

  from dials.algorithms.spot_prediction import PixelLabeller
  labeller = PixelLabeller(beam, detector)

  print 1
  labels = labeller.label(A, 0)
  print 2

  H, K, L = labels.as_vec3_double().parts()
  H.reshape(flex.grid(height, width))
  K.reshape(flex.grid(height, width))
  L.reshape(flex.grid(height, width))

  from matplotlib import pylab
  pylab.imshow(L.as_numpy_array())
  pylab.show()
