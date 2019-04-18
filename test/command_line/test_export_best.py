from __future__ import absolute_import, division, print_function

import procrunner


def test_export_best(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.import",
            "template="
            + dials_data("centroid_test_data").join("centroid_####.cbf").strpath,
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result["exitcode"] and not result["stderr"]
    result = procrunner.run(
        ["dials.find_spots", "imported_experiments.json"],
        working_directory=tmpdir.strpath,
    )
    assert not result["exitcode"] and not result["stderr"]
    result = procrunner.run(
        [
            "dials.index",
            "imported_experiments.json",
            "strong.mpack",
            "space_group=P422",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result["exitcode"] and not result["stderr"]
    result = procrunner.run(
        [
            "dials.integrate",
            "indexed_experiments.json",
            "indexed.mpack",
            "prediction.padding=0",
            "sigma_m_algorithm=basic",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result["exitcode"] and not result["stderr"]
    result = procrunner.run(
        [
            "dials.export",
            "integrated_experiments.json",
            "integrated.mpack",
            "format=best",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result["exitcode"] and not result["stderr"]

    assert tmpdir.join("best.dat").check()
    assert tmpdir.join("best.hkl").check()
    assert tmpdir.join("best.par").check()

    with tmpdir.join("best.dat").open("r") as f:
        lines = "".join(f.readlines()[:10])
    assert (
        lines
        == """\
  183.7743       0.77       1.60
   63.4130       1.57       1.80
   38.3180       1.87       1.71
   27.4540       1.84       1.55
   21.3900       1.89       1.51
   17.5206       1.89       1.52
   14.8370       1.89       1.45
   12.8665       1.90       1.45
   11.3584       1.89       1.42
   10.1669       1.87       1.46
"""
    )

    with tmpdir.join("best.hkl").open("r") as f:
        lines = "".join(f.readlines()[:10])
    assert (
        lines
        == """\
 -20   27   -8      15.92      17.13
 -20   27   -7      64.30      18.73
 -20   27   -6       2.79      16.43
 -20   27   -5      -0.50      16.63
 -20   28  -10      24.24      16.00
 -20   28   -9      46.90      17.03
 -20   28   -7      44.90      18.04
 -20   28   -6      28.86      16.41
 -20   28   -4     -13.18      16.65
 -20   28   -2       7.26      16.54
"""
    )

    lines = tmpdir.join("best.par").read()
    assert (
        lines
        == """\
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
SEPARATION      0.546  0.546
BEAM            219.865  212.610
# end of parameter file for BEST
"""
    )
