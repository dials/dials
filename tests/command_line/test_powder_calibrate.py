from os import path

import pytest

from dials.command_line.powder_calibrate_widget import (
    Geometry,
    PowderCalibrator,
    parse_args,
)


# @pytest.mark.parametrize("eyeball, eyeballed_geom", [(False, "eyeballed.geom")])
def test_pyFAI_from_eyeballed(dials_data):
    test_starting_file = "/home/fyi77748/Data/Al_standard/eyeballed.expt"
    test_args = [
        test_starting_file,
        "standard=Al",
        "eyeball=False",
        "calibrated_geom=/home/fyi77748/Data/Al_standard/test_calibrated.expt",
    ]
    test_expt, test_user = parse_args(args=test_args)
    test_calibrator = PowderCalibrator(expt_params=test_expt, user_args=test_user)
    test_calibrator.calibrate_with_calibrant(verbose=False)
    calibrated_geom = test_calibrator.geometry

    saved_calibration_file = "/home/fyi77748/Data/Al_standard/calibrated.expt"
    saved_args = [saved_calibration_file, "standard=Al", "eyeball=False"]
    saved_expt, _ = parse_args(args=saved_args)
    saved_geom = Geometry(expt_params=saved_expt)

    print((calibrated_geom.param, saved_geom.param))
    assert all(
        [
            pytest.approx(a, 1e-1) == b
            for a, b in zip(calibrated_geom.param, saved_geom.param)
        ]
    )


def test_save_geom_to_expt(dials_data):
    imported_file = "/home/fyi77748/Data/Al_standard/imported.expt"
    imported_args = [imported_file, "standard=Al", "eyeball=False"]
    imported_expt, _ = parse_args(args=imported_args)
    imported_geom = Geometry(expt_params=imported_expt)

    imported_geom.save_to_expt(output="/home/fyi77748/Data/Al_standard/test_save.expt")
    assert path.exists("/home/fyi77748/Data/Al_standard/test_save.expt")

    test_save_file = "/home/fyi77748/Data/Al_standard/test_save.expt"
    test_save_args = [test_save_file, "standard=Al", "eyeball=False"]
    test_save_expt, _ = parse_args(args=test_save_args)
    read_from_saved_geom = Geometry(expt_params=test_save_expt)

    assert imported_geom.param == pytest.approx(read_from_saved_geom.param)
