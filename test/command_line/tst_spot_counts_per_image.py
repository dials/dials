from __future__ import absolute_import, division

def run():
  import os
  import libtbx.load_env
  from libtbx import easy_run
  from glob import glob
  try:
    dials_regression = libtbx.env.dist_path('dials_regression')
  except KeyError:
    print 'FAIL: dials_regression not configured'
    exit(0)

  path = os.path.join(dials_regression, "centroid_test_data")

  # import the data
  cmd = "dials.import %s output.datablock=datablock.json" % ' '.join(glob(os.path.join(path,"*.cbf")))
  easy_run.fully_buffered(cmd).raise_if_errors()
  assert os.path.exists("datablock.json")

  # find the spots
  cmd = "dials.find_spots datablock.json min_spot_size=3"
  easy_run.fully_buffered(cmd).raise_if_errors()
  assert os.path.exists("strong.pickle")

  cmd = "dials.spot_counts_per_image datablock.json strong.pickle plot=spot_counts.png"
  print cmd
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  assert os.path.exists("spot_counts.png"), result.show_stdout()

  assert (
    "| image | #spots | #spots_no_ice | total_intensity |" + \
    " d_min | d_min (distl method 1) | d_min (distl method 2) |"
    in result.stdout_lines), result.stdout_lines


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    run()
