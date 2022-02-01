import pytest

from dials.command_line.powder_calibrate_widget import DialsParams, PowderCalibrator


@pytest.mark.parametrize("eyeball, eyeballed_geom", [(False, "eyeballed.geom")])
def test_powder_calibrate(dials_data, eyeball, eyeballed_geom):
    Al_expt = dials_data("Al_calibrate").join("eyeballed.geom")
    Al_data = DialsParams(expt_file=Al_expt)
    calibrator = PowderCalibrator(Al_data, "Al")
    calibrator.calibrate_with_calibrant()
    calibrated_geom = calibrator.geometry
    print(calibrated_geom)
    assert calibrated_geom
