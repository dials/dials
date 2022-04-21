from __future__ import annotations

import procrunner


def test_export_best(dials_data, tmp_path):
    result = procrunner.run(
        [
            "dials.import",
            "template="
            + str(dials_data("centroid_test_data", pathlib=True) / "centroid_####.cbf"),
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    result = procrunner.run(
        ["dials.find_spots", "imported.expt", "nproc=1"], working_directory=tmp_path
    )
    assert not result.returncode and not result.stderr
    result = procrunner.run(
        ["dials.index", "imported.expt", "strong.refl", "space_group=P422"],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    result = procrunner.run(
        [
            "dials.integrate",
            "nproc=1",
            "indexed.expt",
            "indexed.refl",
            "prediction.padding=0",
            "sigma_m_algorithm=basic",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    result = procrunner.run(
        ["dials.export_best", "integrated.expt", "integrated.refl"],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr

    assert tmp_path.joinpath("best.dat").is_file()
    assert tmp_path.joinpath("best.hkl").is_file()
    assert tmp_path.joinpath("best.par").is_file()

    with tmp_path.joinpath("best.dat").open("r") as f:
        lines = "".join(f.readlines()[:10])
    assert (
        lines
        == """\
  181.8877       0.77       1.60
   63.1895       1.59       1.81
   38.2372       1.87       1.71
   27.4131       1.84       1.55
   21.3655       1.89       1.51
   17.5043       1.88       1.49
   14.8254       1.89       1.45
   12.8580       1.91       1.45
   11.3518       1.89       1.42
   10.1617       1.87       1.41
"""
    )

    with tmp_path.joinpath("best.hkl").open("r") as f:
        lines = "".join(f.readlines()[:10])
    assert (
        lines
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

    lines = tmp_path.joinpath("best.par").read_text()
    assert (
        lines
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
