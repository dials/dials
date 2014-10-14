from __future__ import division

def run():
  import os
  import libtbx.load_env
  from libtbx.test_utils import show_diff
  try:
    dials_regression = libtbx.env.dist_path('dials_regression')
  except KeyError, e:
    print 'FAIL: dials_regression not configured'
    exit(0)

  path = os.path.join(dials_regression, "centroid_test_data")

  # import the data
  from libtbx import easy_run
  cmd = "dials.import %s/*.cbf output=datablock.json" %path
  easy_run.fully_buffered(cmd).raise_if_errors()
  assert os.path.exists("datablock.json")

  # find the spots
  from libtbx import easy_run
  cmd = "dials.find_spots datablock.json min_spot_size=3"
  easy_run.fully_buffered(cmd).raise_if_errors()
  assert os.path.exists("strong.pickle")

  from libtbx import easy_run
  cmd = "dials.spot_counts_per_image datablock.json strong.pickle plot=spot_counts.png"
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  assert os.path.exists("spot_counts.png")

  assert not show_diff("\n".join(result.stdout_lines), """\
The following parameters have been modified:

output {
  plot = spot_counts.png
}
input {
  datablock = datablock.json
  reflections = strong.pickle
}

Per-image analysis:
----------------
|image | #spots|
----------------
|1     | 90    |
|2     | 100   |
|3     | 67    |
|4     | 49    |
|5     | 54    |
|6     | 62    |
|7     | 68    |
|8     | 83    |
|9     | 81    |
----------------\
""")


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    run()
    print "OK"
