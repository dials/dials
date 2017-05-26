from __future__ import absolute_import, division

def run():
  import os
  import libtbx.load_env
  from libtbx import easy_run
  from libtbx.test_utils import show_diff
  try:
    dials_regression = libtbx.env.dist_path('dials_regression')
  except KeyError, e:
    print 'FAIL: dials_regression not configured'
    exit(0)

  path = os.path.join(
    dials_regression, "centroid_test_data/centroid_####.cbf")

  print os.getcwd()

  cmd = "dials.import template=%s" %path
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  cmd = "dials.find_spots datablock.json"
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  cmd = "dials.index datablock.json strong.pickle space_group=P422"
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  cmd = "dials.integrate experiments.json indexed.pickle prediction.padding=0 sigma_m_algorithm=basic"
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  cmd = "dials.export integrated_experiments.json integrated.pickle format=best"
  result = easy_run.fully_buffered(cmd).raise_if_errors()

  assert os.path.exists("best.dat")
  with open("best.dat", "rb") as f:
    lines = ''.join(f.readlines()[:10])
    assert not show_diff(lines, """\
  191.5469       0.00       0.06
   63.8495       1.97       1.39
   38.3104       1.96       1.35
   27.3653       1.59       1.41
   21.2848       1.52       1.41
   17.4156       1.84       2.86
   14.7371       1.86       1.50
   12.7729       1.88       1.47
   11.2710       1.91       2.04
   10.0854       1.86       1.39
""")

  assert os.path.exists("best.hkl")
  with open("best.hkl", "rb") as f:
    lines = ''.join(f.readlines()[:10])
    assert not show_diff(lines, """\
 -20   27   -8      22.61      15.76
 -20   27   -7      69.46      17.54
 -20   27   -6       0.55      15.56
 -20   27   -5      -2.59      15.03
 -20   28  -10      19.99      15.13
 -20   28   -9      38.54      15.95
 -20   28   -7      27.83      16.43
 -20   28   -6      14.81      15.25
 -20   28   -4     -10.41      15.47
 -20   28   -2       6.65      15.26
""")

  assert os.path.exists("best.par")
  with open("best.par", "rb") as f:
    lines = f.read()
    assert not show_diff(lines, """\
# parameter file for BEST
TITLE          From DIALS
DETECTOR       PILA
SITE           Not set
DIAMETER       434.64
PIXEL          0.172
ROTAXIS        -0.01 0.00 1.00 FAST
POLAXIS        0.00 1.00 0.00
GAIN               1.00
CMOSAIC            0.54
PHISTART           0.00
PHIWIDTH           0.20
DISTANCE         191.01
WAVELENGTH      0.97950
POLARISATION    0.99900
SYMMETRY       P422
UB             -0.012247 -0.020072  0.003156
               -0.004952 -0.000405 -0.024640
                0.019676 -0.012595 -0.004237
CELL              42.20    42.20    39.68  90.00  90.00  90.00
RASTER           13  13   7   3   4
SEPARATION      2.960  2.960
BEAM            219.865  212.610
# end of parameter file for BEST
""")
  return


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    run()
    print "OK"
