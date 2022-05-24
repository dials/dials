from __future__ import annotations

import os
from pathlib import Path

import procrunner
import pytest

from dxtbx.serialize import load
from scitbx import matrix

from dials.array_family import flex

base_path = Path(
    os.environ.get(
        "TEST_IMPORT_MTZ", "/dls/mx-scratch/jbe/test_import_mtz_for_multiplex20220520"
    )
)


@pytest.mark.parametrize(
    "pipe, section",
    [
        ("dials", "first"),
        ("dials", "last"),
        ("3dii", "first"),
        ("3dii", "last"),
    ],
)
def test_import_mtz_on_xia2_processing(tmp_path, pipe, section):

    mtz = base_path / f"{pipe}_{section}30.mtz"

    integrated_expt = base_path / f"dials_{section}30_integrated.expt"
    integrated_refl = base_path / f"dials_{section}30_integrated.refl"

    cmd = [
        "dials.import_mtz",
        str(mtz),
        f"output.reflections={section}30.refl",
        f"output.experiments={section}30.expt",
    ]
    result = procrunner.run(cmd, working_directory=tmp_path)
    assert not result.returncode and not result.stderr
    assert (tmp_path / f"{section}30.refl").is_file()
    assert (tmp_path / f"{section}30.expt").is_file()

    imported_expt = load.experiment_list(
        str(tmp_path / f"{section}30.expt"),
        check_format=False,
    )[0]
    imported_refl = flex.reflection_table.from_file(str(tmp_path / f"{section}30.refl"))

    expt_1 = load.experiment_list(str(integrated_expt), check_format=False)[0]
    refl_1 = flex.reflection_table.from_file(str(integrated_refl))

    # Check beam properties. Can't check direction, as MTZ does not store that
    # (?), so it is forced to be along -Z
    assert expt_1.beam.get_wavelength() == pytest.approx(
        imported_expt.beam.get_wavelength()
    )

    # Check detector properties. Can't check precise orientation as MTZ lacks
    # the metadata, so we assume fast along X and slow along -Y.
    assert expt_1.detector[0].get_distance() == pytest.approx(
        imported_expt.detector[0].get_distance(), abs=1.3
    )
    assert expt_1.detector[0].get_pixel_size()[0] == pytest.approx(
        imported_expt.detector[0].get_pixel_size()[0], abs=1e-3
    )
    assert expt_1.detector[0].get_pixel_size()[1] == pytest.approx(
        imported_expt.detector[0].get_pixel_size()[1], abs=1e-3
    )
    assert expt_1.detector[0].get_beam_centre_px(expt_1.beam.get_s0()) == pytest.approx(
        imported_expt.detector[0].get_beam_centre_px(imported_expt.beam.get_s0()),
        abs=4.5,
    )

    # Check scan properties
    assert expt_1.scan.get_image_range() == imported_expt.scan.get_image_range()
    assert expt_1.scan.get_oscillation() == imported_expt.scan.get_oscillation()
    # now for the crystal, as we are trying to put everything into our DIALS
    # geometry, both the U and B matrices should be the same
    assert expt_1.crystal.get_B() == pytest.approx(
        imported_expt.crystal.get_B(), abs=2e-3
    )
    # FIXME
    # assert expt_1.crystal.get_U() == pytest.approx(
    #    imported_expt.crystal.get_U(), abs=1e-2
    # )

    # Now test phi etc
    phi = imported_refl["xyzobs.px.value"].parts()[2]
    if section == "first":
        assert max(phi) < 30 and min(phi) > 0
    else:
        assert max(phi) < 90 and min(phi) > 60

    # choose a few random reflections to compare
    check_last = [(8, 18, 11), (-16, 20, 23), (8, 20, 10)]
    check_first = [(5, 3, 22), (-1, 5, 3), (-3, -37, -7)]
    if pipe == "dials":
        if section == "last":
            miller_index_to_check = check_last
        else:
            miller_index_to_check = check_first
        for idx in miller_index_to_check:
            r1 = refl_1.select(refl_1["miller_index"] == idx)
            r1 = r1.select(r1.get_flags(r1.flags.integrated_prf))
            r2 = imported_refl.select(imported_refl["miller_index"] == idx)
            # for these examples, there is only one reflection
            assert r1[0]["intensity.prf.value"] == pytest.approx(
                r2[0]["intensity.prf.value"], abs=1e-2
            )
            # Can't reasonably check absolute directions, but can check
            # scattering angle
            tt1 = matrix.col(expt_1.beam.get_s0()).angle(
                matrix.col(r1[0]["s1"]), deg=True
            )
            tt2 = matrix.col(imported_expt.beam.get_s0()).angle(
                matrix.col(r2[0]["s1"]), deg=True
            )
            assert tt1 == pytest.approx(tt2, abs=0.4)  # abs is surprisingly high here!

            # note in dials, we export the observed x an y, and the calculated z.
            assert r1[0]["xyzcal.px"][2] == pytest.approx(
                r2[0]["xyzcal.px"][2], abs=2e-2
            )
            assert r1[0]["xyzobs.px.value"][0] == pytest.approx(
                r2[0]["xyzobs.px.value"][0], abs=2e-2
            )
            assert r1[0]["xyzobs.px.value"][1] == pytest.approx(
                r2[0]["xyzobs.px.value"][1], abs=2e-2
            )


@pytest.mark.parametrize("pipe", ["dials", "3dii"])
def test_multiplex_on_imported_mtz(tmp_path, pipe):

    for section in ["first", "last"]:
        cmd = [
            "dials.import_mtz",
            str(base_path / f"{pipe}_{section}30.mtz"),
            f"output.reflections=imported_{section}.refl",
            f"output.experiments=imported_{section}.expt",
        ]
        result = procrunner.run(cmd, working_directory=tmp_path)
        assert not result.returncode and not result.stderr

    cmd = ["xia2.multiplex", "unit_cell.refine=None"]
    for section in ["first", "last"]:
        cmd.extend(
            [
                str(tmp_path / f"imported_{section}.expt"),
                str(tmp_path / f"imported_{section}.refl"),
            ]
        )

    result = procrunner.run(cmd, working_directory=tmp_path)
    assert not result.returncode and not result.stderr


"""


base_path = os.environ.get("TEST_IMPORT_MTZ", "/dls/mx-scratch/jbe/test_import_mtz/")
dials_expt = base_path + "dials_integrated.expt"
dials_refl = base_path + "dials_integrated.refl"


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

### how about an old mosflm file?"""
