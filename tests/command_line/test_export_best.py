from __future__ import annotations

import shutil
import subprocess


def test_export_best(dials_data, tmp_path):
    # make sure the raw image data is available as it is required for export_best
    _ = dials_data("insulin")
    integrated_expt = dials_data("insulin_processed") / "integrated.expt"
    integrated_refl = dials_data("insulin_processed") / "integrated.refl"
    result = subprocess.run(
        [shutil.which("dials.export_best"), integrated_expt, integrated_refl],
        cwd=tmp_path,
        capture_output=True,
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
  270.0172       8.91      10.06
   92.2098      15.00       9.17
   55.5986      17.71       7.67
   39.7977      28.05      15.20
   30.9906      71.74      28.53
   25.3755     146.49      49.59
   21.4832     200.23      52.15
   18.6264     216.27      51.17
   16.4404     226.93      50.73
   14.7137     234.22      49.41
"""
    )

    with tmp_path.joinpath("best.hkl").open("r") as f:
        lines = "".join(f.readlines()[:10])
    assert (
        lines
        == """\
 -43  -17   28    1971.02     101.36
 -43  -16   27    1489.20     141.01
 -43  -15   26     589.66     173.22
 -43  -15   28    1161.59     142.63
 -43  -14   23     -19.91     100.28
 -43  -14   25     -71.86     140.74
 -43  -14   27    -345.97     174.45
 -43  -14   29      -5.58     141.07
 -43  -13   22    -191.52     141.15
 -43  -13   24    -172.19     172.15
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
DIAMETER       188.01
PIXEL          0.0816
ROTAXIS        0.00 0.00 1.00 FAST
POLAXIS        0.00 1.00 0.00
GAIN               1.00
CMOSAIC            0.35
PHISTART           0.00
PHIWIDTH           1.00
DISTANCE         158.71
WAVELENGTH      0.97900
POLARISATION    0.99900
SYMMETRY       I23
UB              0.001789 -0.007488 -0.010232
                0.001789  0.010379 -0.007283
                0.012553 -0.000412  0.002496
CELL              78.09    78.09    78.09  90.00  90.00  90.00
RASTER           7 7 5 3 3
SEPARATION      0.448  0.448
BEAM             94.415   94.513
# end of parameter file for BEST
"""
    )
