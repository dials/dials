import pytest

from dials.command_line.powder_calibrate_widget import DialsParams, PowderCalibrator


@pytest.mark.parametrize("eyeball, eyeballed_geom", [(False, "eyeballed.geom")])
def test_powder_calibrate(dials_data, eyeball, eyeballed_geom):
    # Al_expt = dials_data("Al_calibrate").join("eyeballed.geom")
    Al_expt = "/home/elena/Work/Diamond/Data/Al_data/pyFAI/dials_data/imported.expt"
    eyeball_file = (
        "/home/elena/Work/Diamond/Data/Al_data/pyFAI/dials_data/eyeballed.geom"
    )
    Al_data = DialsParams(expt_file=Al_expt)
    calibrator = PowderCalibrator(
        Al_data, standard="Al", eyeball=False, start_geom_file=eyeball_file
    )
    calibrator.calibrate_with_calibrant(verbose=False)
    calibrated_geom = calibrator.geometry

    assert calibrated_geom
