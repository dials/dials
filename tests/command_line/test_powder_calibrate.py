from os import path
from unittest.mock import patch

import pytest

from dials.command_line.powder_calibrate_widget import (
    Geometry,
    PowderCalibrator,
    parse_args,
)


@patch("matplotlib.pyplot.show")
@patch("matplolib.widget.Slider.on_changed")
def test_calibrate_coarse(dials_data, tmpdir, eyeball, starting_geometry):
    aluminium_powder = dials_data("aluminium_standard", pathlib=True)
    starting_geom = aluminium_powder / starting_geometry

    test_calibrator = PowderCalibrator(
        starting_geom,
        standard="Al",
        eyeball=eyeball,
        calibrated_geom=str(tmpdir / "test_calibrated.expt"),
    )

    with patch.object(
        test_calibrator,
        "val",
    ):
        test_calibrator.calibrate_with_calibrant(plots=False)
    calibrated_geom = test_calibrator.geometry

    expected_calibrated_file = aluminium_powder / "calibrated.expt"
    expected_expt, _ = parse_args(
        args=[str(expected_calibrated_file), "standard=Al", "eyeball=False"]
    )
    expected_geom = Geometry(expt_params=expected_expt)

    assert all(
        [
            pytest.approx(a, 1e-1) == b
            for a, b in zip(calibrated_geom.param, expected_geom.param)
        ]
    )


def test_save_geom_to_expt(dials_data, tmpdir):
    aluminium_powder = dials_data("aluminium_standard", pathlib=True)
    imported = aluminium_powder / "imported.expt"
    imported_expt, _ = parse_args(args=[str(imported), "standard=Al", "eyeball=False"])
    imported_geom = Geometry(expt_params=imported_expt)
    outfile = str(tmpdir) + "/test_save.expt"
    imported_geom.save_to_expt(output=outfile)
    assert path.exists(outfile)

    test_save_expt, _ = parse_args(args=[outfile, "standard=Al", "eyeball=False"])
    read_from_saved_geom = Geometry(expt_params=test_save_expt)

    assert imported_geom.param == pytest.approx(read_from_saved_geom.param)
