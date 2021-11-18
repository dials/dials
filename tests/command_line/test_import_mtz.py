import os

import numpy as np
import procrunner
import pytest

from dxtbx.serialize import load

from dials.array_family import flex

base_path = os.environ.get("TEST_IMPORT_MTZ", "/dls/mx-scratch/jbe/test_import_mtz/")
dials_expt = base_path + "dials_integrated.expt"
dials_refl = base_path + "dials_integrated.refl"


def check_dials_refls(integrated_refls, imported_refls):
    # now check some random reflections
    miller_index_to_check = [(9, 7, 3), (21, 13, 4), (8, -26, -12)]
    for idx in miller_index_to_check:
        r1 = integrated_refls.select(integrated_refls["miller_index"] == idx)
        r1 = r1.select(r1.get_flags(r1.flags.integrated_prf))
        r2 = imported_refls.select(imported_refls["miller_index"] == idx)
        # for these examples, there is only one reflection
        assert r1[0]["intensity.prf.value"] == pytest.approx(
            r2[0]["intensity.prf.value"], abs=1e-2
        )
        assert r1[0]["s1"] == pytest.approx(r2[0]["s1"], abs=1e-2)
        # note in dials, we export the observed x an y, and the calculated z.
        assert r1[0]["xyzcal.px"][2] == pytest.approx(r2[0]["xyzcal.px"][2], abs=1e-2)
        assert r1[0]["xyzobs.px.value"][0] == pytest.approx(
            r2[0]["xyzobs.px.value"][0], abs=1e-2
        )
        assert r1[0]["xyzobs.px.value"][1] == pytest.approx(
            r2[0]["xyzobs.px.value"][1], abs=1e-2
        )


def test_import_dials_integrated_mtz_with_template(dials_data, tmp_path):
    integrated_mtz = base_path + "dials_INTEGRATE.mtz"
    image_template = os.path.join(
        dials_data("x4wide").strpath, "X4_wide_M1S4_2_####.cbf"
    )

    result = procrunner.run(
        [
            "dials.import_mtz",
            integrated_mtz,
            f"images.template={image_template}",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "imported_mtz.expt").is_file()
    assert (tmp_path / "imported_mtz.refl").is_file()

    imported_expts = load.experiment_list(
        tmp_path / "imported_mtz.expt", check_format=False
    )
    final_expts = load.experiment_list(dials_expt, check_format=False)

    assert imported_expts[0].goniometer == final_expts[0].goniometer
    assert imported_expts[0].beam.get_wavelength() == pytest.approx(
        final_expts[0].beam.get_wavelength()
    )
    assert imported_expts[0].beam.get_unit_s0() == pytest.approx(
        final_expts[0].beam.get_unit_s0(), abs=1e-2
    )  # beam is refined in dials
    # note don't compare detectors, several components refined
    assert (
        imported_expts[0].scan.get_image_range()
        == final_expts[0].scan.get_image_range()
    )
    assert (
        imported_expts[0].scan.get_oscillation()
        == final_expts[0].scan.get_oscillation()
    )
    # now for the crystal, as we are trying to put everything into our DIALS
    # geometry, both the U and B matrices should be the same
    assert imported_expts[0].crystal.get_B() == pytest.approx(
        final_expts[0].crystal.get_B(), abs=1e-3
    )
    assert imported_expts[0].crystal.get_U() == pytest.approx(
        final_expts[0].crystal.get_U(), abs=1e-2
    )

    integrated_refls = flex.reflection_table.from_file(dials_refl)
    imported_refls = flex.reflection_table.from_file(tmp_path / "imported_mtz.refl")
    check_dials_refls(integrated_refls, imported_refls)


def test_import_dials_integrated_mtz_without_template(tmp_path):
    integrated_mtz = base_path + "dials_INTEGRATE.mtz"

    result = procrunner.run(
        [
            "dials.import_mtz",
            integrated_mtz,
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "imported_mtz.expt").is_file()
    assert (tmp_path / "imported_mtz.refl").is_file()

    imported_expts = load.experiment_list(
        tmp_path / "imported_mtz.expt", check_format=False
    )
    final_expts = load.experiment_list(dials_expt, check_format=False)

    assert imported_expts[0].goniometer == final_expts[0].goniometer
    assert imported_expts[0].beam.get_wavelength() == pytest.approx(
        final_expts[0].beam.get_wavelength()
    )
    assert imported_expts[0].beam.get_unit_s0() == pytest.approx(
        final_expts[0].beam.get_unit_s0(), abs=1e-2
    )  # beam is refined in dials
    assert not imported_expts[0].imageset
    # note don't compare detectors, several components refined
    # for the scan, the only item not set correctly is the exposure time
    assert (
        imported_expts[0].scan.get_image_range()
        == final_expts[0].scan.get_image_range()
    )
    assert (
        imported_expts[0].scan.get_oscillation()
        == final_expts[0].scan.get_oscillation()
    )

    # now for the crystal, as we are trying to put everything into our DIALS
    # geometry, both the U and B matrices should be the same
    assert imported_expts[0].crystal.get_B() == pytest.approx(
        final_expts[0].crystal.get_B(), abs=1e-3
    )
    assert imported_expts[0].crystal.get_U() == pytest.approx(
        final_expts[0].crystal.get_U(), abs=1e-2
    )

    integrated_refls = flex.reflection_table.from_file(dials_refl)
    imported_refls = flex.reflection_table.from_file(tmp_path / "imported_mtz.refl")
    check_dials_refls(integrated_refls, imported_refls)


def test_import_3dii_integrated_mtz_with_template(dials_data, tmp_path):
    integrated_mtz = base_path + "3dii_INTEGRATE.mtz"
    image_template = os.path.join(
        dials_data("x4wide").strpath, "X4_wide_M1S4_2_####.cbf"
    )

    result = procrunner.run(
        [
            "dials.import_mtz",
            integrated_mtz,
            f"images.template={image_template}",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "imported_mtz.expt").is_file()
    assert (tmp_path / "imported_mtz.refl").is_file()

    imported_expts = load.experiment_list(
        tmp_path / "imported_mtz.expt", check_format=False
    )
    final_expts = load.experiment_list(dials_expt, check_format=False)

    assert imported_expts[0].goniometer == final_expts[0].goniometer
    assert imported_expts[0].beam.get_wavelength() == pytest.approx(
        final_expts[0].beam.get_wavelength()
    )
    assert imported_expts[0].beam.get_unit_s0() == pytest.approx(
        final_expts[0].beam.get_unit_s0(), abs=1e-2
    )  # beam is refined in dials
    # assert imported_expts[0].imageset == final_expts[0].imageset
    # note don't compare detectors, several components refined
    assert (
        imported_expts[0].scan.get_image_range()
        == final_expts[0].scan.get_image_range()
    )
    assert (
        imported_expts[0].scan.get_oscillation()
        == final_expts[0].scan.get_oscillation()
    )
    # now for the crystal, as we are trying to put everything into our DIALS
    # geometry, both the U and B matrices should be the same
    assert imported_expts[0].crystal.get_B() == pytest.approx(
        final_expts[0].crystal.get_B(), abs=1e-3
    )
    # for the U matrix, the crystal parameters in this example are sufficiently
    # different that some elements are not so similar, but on the whole they
    # should be similar enough, so multiply U1T * U2 and see if diagonal is
    # similar to the identity matrix.
    U1 = np.array(imported_expts[0].crystal.get_U()).reshape(3, 3)
    U2 = np.array(final_expts[0].crystal.get_U()).reshape(3, 3)
    multiplication = U1.T @ U2
    assert all(multiplication[i, i] > 0.96 for i in range(3))


def test_import_3dii_integrated_mtz_without_template(tmp_path):
    integrated_mtz = base_path + "3dii_INTEGRATE.mtz"

    result = procrunner.run(
        [
            "dials.import_mtz",
            integrated_mtz,
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "imported_mtz.expt").is_file()
    assert (tmp_path / "imported_mtz.refl").is_file()

    imported_expts = load.experiment_list(
        tmp_path / "imported_mtz.expt", check_format=False
    )
    final_expts = load.experiment_list(dials_expt, check_format=False)

    assert imported_expts[0].goniometer == final_expts[0].goniometer
    assert imported_expts[0].beam.get_wavelength() == pytest.approx(
        final_expts[0].beam.get_wavelength()
    )
    assert imported_expts[0].beam.get_unit_s0() == pytest.approx(
        final_expts[0].beam.get_unit_s0(), abs=1e-2
    )  # beam is refined in dials
    assert not imported_expts[0].imageset
    # note don't compare detectors, several components refined
    # for the scan, the only item not set correctly is the exposure time
    assert (
        imported_expts[0].scan.get_image_range()
        == final_expts[0].scan.get_image_range()
    )
    assert (
        imported_expts[0].scan.get_oscillation()
        == final_expts[0].scan.get_oscillation()
    )
    # now for the crystal, as we are trying to put everything into our DIALS
    # geometry, both the U and B matrices should be the same
    assert imported_expts[0].crystal.get_B() == pytest.approx(
        final_expts[0].crystal.get_B(), abs=1e-3
    )
    # for the U matrix, the crystal parameters in this example are sufficiently
    # different that some elements are not so similar, but on the whole they
    # should be similar enough, so multiply U1T * U2 and see if diagonal is
    # similar to the identity matrix.
    U1 = np.array(imported_expts[0].crystal.get_U()).reshape(3, 3)
    U2 = np.array(final_expts[0].crystal.get_U()).reshape(3, 3)
    multiplication = U1.T @ U2
    assert all(multiplication[i, i] > 0.96 for i in range(3))


### test a suitable error is raised if a subset of images have been used and
### there is a mismatch with the scan range.

### separate test for data processed purely with xds (not 3dii)?

### how about an old mosflm file?
