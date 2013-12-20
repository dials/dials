from __future__ import division

import os
import cPickle as pickle
from libtbx.test_utils import open_tmp_file, show_diff
from libtbx import easy_run

from dials.model.data import ReflectionList # import dependency

def exercise_spots_xds():
  xparm_txt = """\
 XPARM.XDS
     1       -5.0000    0.2000  1.000000  0.000000  0.000000
       1.252370      -0.000017       0.001994       0.798484
   143     93.7855     93.7855    129.9110  90.000  90.000 120.000
      36.211739     -85.308029     -14.386198
     -34.227608      49.189667     -72.142654
     117.029358      52.951710     -19.419268
         1      2463      2527    0.172000    0.172000
    1226.947998    1222.105469     180.111465
       1.000000       0.000000       0.000000
       0.000000       1.000000       0.000000
       0.000000       0.000000       1.000000
         1         1      2463         1      2527
    0.00    0.00    0.00  1.00000  0.00000  0.00000  0.00000  1.00000  0.00000
"""

  f = open_tmp_file(suffix="_XPARM.XDS", mode="wb")
  f.write(xparm_txt)
  f.close()
  xparm_file = f.name

  txt = """\
 2411.40 1000.70 25.00 16384. 0 0 0
 1328.60 2170.40 20.57 7326. 0 0 0
 177.56 2191.30 24.94 6779. 0 0 0
 1231.34 1220.04 24.99 1952. 0 0 0
 1227.07 1230.56 24.81 349. 0 0 0
 1341.63 1243.25 5.64 321. 2 -2 11
 1125.23 1197.72 12.14 231. -1 2 -10
 1317.52 1171.59 19.28 120. 6 -4 6
 1260.25 1300.55 13.67 116. -4 2 6
 1090.27 1199.47 41.49 114. -2 3 -13
"""

  f = open_tmp_file(suffix="_SPOTS.XDS", mode="wb")
  f.write(txt)
  f.close()

  output_pickle = "%s.pickle" %f.name[:-4]
  args = ["dials.extract_spot_xds", f.name, xparm_file,
          "-o '%s'" %output_pickle]
  result = easy_run.fully_buffered(command=" ".join(args)).raise_if_errors()
  assert os.path.exists(output_pickle)
  with open(output_pickle, "rb") as f:
    reflections = pickle.load(f)
    assert len(reflections) == 5

  # now test we can export it
  args = ["dials.export_spot_xds", output_pickle]
  result = easy_run.fully_buffered(command=" ".join(args)).raise_if_errors()
  assert os.path.exists("SPOT.XDS")
  with open("SPOT.XDS", "rb") as f:
    txt = f.read()
    assert not show_diff(
      "\n".join([l.rstrip() for l in txt.split("\n")]), """\
 230.67 213.75 -0.07 321.00  2 -2 11
 193.45 205.92 -0.04 231.00  -1 2 -10
 226.53 201.43 -0.02 120.00  6 -4 6
 216.68 223.61 -0.04 116.00  -4 2 6
 187.44 206.22 0.06 114.00  -2 3 -13
""")


def export_xds():
  from libtbx import easy_run
  import libtbx.load_env

  try:
    dials_regression = libtbx.env.dist_path('dials_regression')
  except KeyError, e:
    print "skipping exercise_export_xds: dials_regression not available"
    return

  args = ["dials.export_xds",
          os.path.join(dials_regression, "centroid_test_data",
                       "crystal.json"),
          os.path.join(dials_regression, "centroid_test_data",
                       "sweep.json"),
          os.path.join(dials_regression, "centroid_test_data",
                       "spot_all_xds.pickle")]
  result = easy_run.fully_buffered(command=" ".join(args)).raise_if_errors()
  assert os.path.exists("XDS.INP")
  assert os.path.exists("XPARM.XDS")
  assert os.path.exists("SPOT.XDS")


def run():
  export_xds()
  exercise_spots_xds()

if __name__ == '__main__':
  run()
  print "OK"
