from __future__ import division
import os
import cPickle as pickle
import libtbx.load_env
from libtbx import easy_run
from libtbx.test_utils import approx_equal
from dials.model.data import ReflectionList # import dependency

def exercise_spotfinder():
  if not libtbx.env.has_module("dials_regression"):
    print "Skipping exercise_spotfinder: dials_regression not present"
    return

  data_dir = libtbx.env.find_in_repositories(
    relative_path="dials_regression/centroid_test_data",
    test=os.path.isdir)
  template = os.path.join(data_dir, "centroid_00*.cbf")
  args = ["dials.spotfinder", template, "-o spotfinder.pickle", "--nproc=1"]
  result = easy_run.fully_buffered(command=" ".join(args)).raise_if_errors()
  assert os.path.exists("spotfinder.pickle")
  with open("spotfinder.pickle", "rb") as f:
    reflections = pickle.load(f)
    assert len(reflections) == 459
    refl = reflections[0]
    assert approx_equal(refl.intensity, 142)
    assert approx_equal(refl.frame_number, 0)
    assert approx_equal(refl.bounding_box, (1258, 1260, 537, 541, 0, 1))
    assert approx_equal(refl.centroid_position,
                        (1258.7957746478874, 539.112676056338, 0.5))

  # now with a resolution filter
  args = ["dials.spotfinder", "d_min=2", "d_max=15",
          template, "-o spotfinder.pickle", "--nproc=1"]
  result = easy_run.fully_buffered(command=" ".join(args)).raise_if_errors()
  assert os.path.exists("spotfinder.pickle")
  with open("spotfinder.pickle", "rb") as f:
    reflections = pickle.load(f)
    assert len(reflections) == 371

  # now with more generous parameters
  args = ["dials.spotfinder", "min_spot_size=3", "max_separation=3",
          template, "-o spotfinder.pickle", "--nproc=1"]
  result = easy_run.fully_buffered(command=" ".join(args)).raise_if_errors()
  assert os.path.exists("spotfinder.pickle")
  with open("spotfinder.pickle", "rb") as f:
    reflections = pickle.load(f)
    assert len(reflections) == 679

  # now with XFEL stills
  data_dir = libtbx.env.find_in_repositories(
    relative_path="dials_regression/spotfinding_test_data",
    test=os.path.isdir)
  template = os.path.join(data_dir, "idx*.cbf")
  args = ["dials.spotfinder", template, "-o spotfinder.pickle", "--nproc=1"]
  result = easy_run.fully_buffered(command=" ".join(args)).raise_if_errors()
  assert os.path.exists("spotfinder.pickle")
  with open("spotfinder.pickle", "rb") as f:
    reflections = pickle.load(f)
    assert len(reflections) == 2989

def run():
  exercise_spotfinder()

if __name__ == '__main__':
  run()
  print 'OK'
