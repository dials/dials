from __future__ import absolute_import, division
import os
import cPickle as pickle
import libtbx.load_env
from libtbx import easy_run
from dials.array_family import flex # import dependency

def exercise_merge_spot_lists():
  if not libtbx.env.has_module("dials_regression"):
    print "Skipping exercise_spotfinder: dials_regression not present"
    return

  data_dir = libtbx.env.find_in_repositories(
    relative_path="dials_regression/centroid_test_data",
    test=os.path.isdir)

  template = os.path.join(data_dir, "centroid_00*.cbf")
  args = ["dials.spotfinder", template, "-o spotfinder.pickle",
    "min_spot_size=0", "max_separation=100"]
  result = easy_run.fully_buffered(command=" ".join(args)).raise_if_errors()
  assert os.path.exists("spotfinder.pickle")

  template = os.path.join(data_dir, "centroid_00*.cbf")
  args = ["dials.spotfinder", template, "-o spotfinder1.pickle",
    "min_spot_size=0", "max_separation=100", 'scan_range="0, 4"']
  result = easy_run.fully_buffered(command=" ".join(args)).raise_if_errors()
  assert os.path.exists("spotfinder1.pickle")

  template = os.path.join(data_dir, "centroid_00*.cbf")
  args = ["dials.spotfinder", template, "-o spotfinder2.pickle",
    "min_spot_size=0", "max_separation=100", 'scan_range="4, 9"']
  result = easy_run.fully_buffered(command=" ".join(args)).raise_if_errors()
  assert os.path.exists("spotfinder2.pickle")

  args = ["dials.merge_spot_lists", 'spotfinder1.pickle', 'spotfinder2.pickle',
    "-o spotfinder3.pickle"]
  result = easy_run.fully_buffered(command=" ".join(args)).raise_if_errors()
  assert os.path.exists("spotfinder3.pickle")

  spot1 = pickle.load(open('spotfinder.pickle'))
  spot2 = pickle.load(open('spotfinder3.pickle'))
  assert(len(spot1) == len(spot2))

  for r1, r2 in zip(spot1, spot2):
    from scitbx import matrix
    c1 = matrix.col(r1.centroid_position)
    c2 = matrix.col(r2.centroid_position)
    assert(abs(c1 - c2) < 1e-7)
    i1 = r1.intensity
    i2 = r2.intensity
    assert(abs(i1 - i2) < 1e-7)

  print 'OK'

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    exercise_merge_spot_lists()
