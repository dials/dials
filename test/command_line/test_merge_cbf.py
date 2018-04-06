from __future__ import absolute_import, division, print_function

import glob
import os

from libtbx import easy_run
from libtbx.test_utils import open_tmp_directory

def test_merge_cbf(dials_regression, tmpdir):
  tmpdir.chdir()
  data_dir = os.path.join(dials_regression, "centroid_test_data")

  g = glob.glob(os.path.join(data_dir, "*.cbf"))
  assert len(g) == 9

  cmd = "dials.merge_cbf %s merge_n_images=3" %(" ".join(g))
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  g = glob.glob(os.path.join(tmpdir.strpath, "sum_*.cbf"))
  assert len(g) == 3

  # test alternate mode of accessing image data
  cmd += " image_prefix=sum2_ get_raw_data_from_imageset=false"
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  g2 = glob.glob(os.path.join(tmpdir.strpath, "sum2_*.cbf"))
  assert len(g2) == 3

  # check summed images are the same in either case
  import filecmp
  for f1, f2 in zip(sorted(g), sorted(g2)):
    assert filecmp.cmp(f1, f2)
