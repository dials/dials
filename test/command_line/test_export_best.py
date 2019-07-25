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
        ["dials.find_spots", "imported.expt"], working_directory=tmpdir.strpath
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
        ["dials.export", "integrated.expt", "integrated.refl", "format=best"],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr

    assert tmpdir.join("best.dat").check()
    assert tmpdir.join("best.hkl").check()
    assert tmpdir.join("best.par").check()

    with tmpdir.join("best.dat").open("r") as f:
        lines = "".join(f.readlines()[:10])
    assert (
        lines
        == """\
  183.7552       0.77       1.60
   63.4094       1.57       1.80
   38.3162       1.87       1.71
   27.4528       1.84       1.55
   21.3892       1.89       1.51
   17.5199       1.89       1.52
   14.8364       1.89       1.45
   12.8661       1.90       1.45
   11.3580       1.89       1.42
   10.1665       1.87       1.46
"""
    )

    with tmpdir.join("best.hkl").open("r") as f:
        lines = "".join(f.readlines()[:10])
    assert (
        lines
        == """\
 -36   -4   17       3.93      13.67
 -36   -4   18     -10.54      13.65
 -36   -3   16      18.34      16.24
 -36   -2   15      -1.53      16.57
 -36   -1   13     -14.19      17.01
 -36    0   11     -22.41      16.27
 -36    0   12      -2.57      14.54
 -36    1    9      25.03      16.88
 -36    1   10      -8.71      15.93
 -35   -4   18      29.57      17.44
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
ROTAXIS        -0.00 0.00 1.00 FAST
POLAXIS        0.00 1.00 0.00
GAIN               1.00
CMOSAIC            0.54
PHISTART           0.00
PHIWIDTH           0.20
DISTANCE         190.99
WAVELENGTH      0.97950
POLARISATION    0.99900
SYMMETRY       P422
UB              0.012248  0.020072  0.003156
                0.005026  0.000358 -0.024626
               -0.019659  0.012598 -0.004329
CELL              42.19    42.19    39.68  90.00  90.00  90.00
RASTER           7 7 5 3 3
SEPARATION      0.545  0.545
BEAM            219.865  212.610
# end of parameter file for BEST
"""
    )
