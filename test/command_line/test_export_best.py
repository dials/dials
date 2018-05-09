from __future__ import absolute_import, division, print_function

import os
import procrunner

def test_export_best(dials_regression, tmpdir):
  tmpdir.chdir()
  path = os.path.join(
    dials_regression, "centroid_test_data", "centroid_####.cbf")

  result = procrunner.run(["dials.import", "template=" + path])
  assert not result['exitcode'] and not result['stderr']
  result = procrunner.run(["dials.find_spots", "datablock.json"])
  assert not result['exitcode'] and not result['stderr']
  result = procrunner.run(["dials.index", "datablock.json", "strong.pickle", "space_group=P422"])
  assert not result['exitcode'] and not result['stderr']
  result = procrunner.run([
      "dials.integrate",
      "experiments.json",
      "indexed.pickle",
      "prediction.padding=0",
      "sigma_m_algorithm=basic",
  ])
  assert not result['exitcode'] and not result['stderr']
  result = procrunner.run(["dials.export", "integrated_experiments.json", "integrated.pickle", "format=best"])
  assert not result['exitcode'] and not result['stderr']

  assert os.path.exists("best.dat")
  assert os.path.exists("best.hkl")
  assert os.path.exists("best.par")

  with open("best.dat", "r") as f:
    lines = ''.join(f.readlines()[:10])
  assert lines == """\
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
"""

  with open("best.hkl", "r") as f:
    lines = ''.join(f.readlines()[:10])
  assert lines == """\
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
"""

  with open("best.par", "r") as f:
    lines = f.read()
  assert lines == """\
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
RASTER           7 7 5 3 3
SEPARATION      0.500  0.500
BEAM            219.865  212.610
# end of parameter file for BEST
"""
