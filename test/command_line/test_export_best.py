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
    assert not result.returncode and not result.stderr
    result = procrunner.run(
        ["dials.find_spots", "imported.expt", "nproc=1"],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    result = procrunner.run(
        ["dials.index", "imported.expt", "strong.refl", "space_group=P422"],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    result = procrunner.run(
        [
            "dials.integrate",
            "indexed.expt",
            "indexed.refl",
            "prediction.padding=0",
            "sigma_m_algorithm=basic",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    result = procrunner.run(
        ["dials.export_best", "integrated.expt", "integrated.refl"],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr

    assert tmpdir.join("best.dat").check()
    assert tmpdir.join("best.hkl").check()
    assert tmpdir.join("best.par").check()

    with tmpdir.join("best.dat").open("r") as f:
        lines = "".join(f.readlines()[:10])
    assert (
        lines.replace("\r", "")
        == """\
  181.8877       0.50       1.60
   63.1895       0.63       1.90
   38.2372       1.15       1.93
   27.4131       1.41       1.75
   21.3655       1.57       1.69
   17.5043       1.62       1.65
   14.8254       1.68       1.59
   12.8580       1.72       1.57
   11.3518       1.72       1.53
   10.1617       1.72       1.51
"""
    )

    with tmpdir.join("best.hkl").open("r") as f:
        lines = "".join(f.readlines()[:10])
    assert (
        lines.replace("\r", "")
        == """\
 -20   27   -8      20.17      20.00
 -20   27   -7      74.13      21.59
 -20   27   -6      22.34      19.57
 -20   27   -5       6.33      19.72
 -20   28  -10      19.77      18.77
 -20   28   -9      50.37      20.28
 -20   28   -7      69.23      21.42
 -20   28   -6      24.42      19.56
 -20   28   -4     -10.35      19.51
 -20   28   -2      47.53      20.49
"""
    )

    lines = tmpdir.join("best.par").read()
    assert (
        lines.replace("\r", "")
        == """\
# parameter file for BEST
TITLE          From DIALS
DETECTOR       PILA
SITE           Not set
DIAMETER       434.64
PIXEL          0.172
ROTAXIS        0.01 0.00 1.00 FAST
POLAXIS        0.00 1.00 0.00
GAIN               1.00
CMOSAIC            0.54
PHISTART           0.00
PHIWIDTH           0.20
DISTANCE         191.09
WAVELENGTH      0.97950
POLARISATION    0.99900
SYMMETRY       P422
UB             -0.012248 -0.020067  0.003152
               -0.005029 -0.000351 -0.024623
                0.019651 -0.012597 -0.004336
CELL              42.20    42.20    39.68  90.00  90.00  90.00
RASTER           7 7 5 3 3
SEPARATION      0.665  0.665
BEAM            219.875  212.612
# end of parameter file for BEST
"""
    )
