from __future__ import absolute_import, division, print_function

from glob import glob
import os
import procrunner

def test(dials_regression, run_in_tmpdir):
  images = glob(os.path.join(dials_regression, "centroid_test_data", "centroid*.cbf"))

  result = procrunner.run_process([
      "dials.find_spots",
      "output.experiments=experiments.json",
      "output.reflections=spotfinder.pickle",
      "output.shoeboxes=True",
  ] + images)
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("experiments.json")
  assert os.path.exists("spotfinder.pickle")

  result = procrunner.run_process([
      "dials.find_hot_pixels",
      "input.experiments=experiments.json",
      "input.reflections=spotfinder.pickle",
      "output.mask=hot_mask.pickle"
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("hot_mask.pickle")
  assert "Found 8 hot pixels" in result['stdout']
