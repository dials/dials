from __future__ import absolute_import, division, print_function

import os

import procrunner
import pytest


def test_multiple_sweep_import_fails_without_allow_parameter(dials_data, tmpdir):
    # Find the image files
    image_files = dials_data("centroid_test_data").listdir("centroid*.cbf", sort=True)
    del image_files[4]  # Delete filename to force two sweeps

    # run without allowing multiple sweeps
    result = procrunner.run(
        ["dials.import", "output.experiments=experiments_multiple_sweeps.expt"]
        + [f.strpath for f in image_files],
        working_directory=tmpdir.strpath,
    )
    assert result.returncode == 1
    assert b"ore than 1 sweep" in result.stderr
    assert not tmpdir.join("experiments_multiple_sweeps.expt").check()


def test_multiple_sweep_import_suceeds_with_allow_parameter(dials_data, tmpdir):
    # Find the image files
    image_files = dials_data("centroid_test_data").listdir("centroid*.cbf", sort=True)
    del image_files[4]  # Delete filename to force two sweeps

    result = procrunner.run(
        [
            "dials.import",
            "output.experiments=experiments_multiple_sweeps.expt",
            "allow_multiple_sweeps=True",
        ]
        + [f.strpath for f in image_files],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("experiments_multiple_sweeps.expt").check(file=1)

    from dxtbx.serialize import load

    experiments = load.experiment_list(
        tmpdir.join("experiments_multiple_sweeps.expt").strpath
    )
    assert len(experiments) == 2


def test_with_mask(dials_data, tmpdir):
    # Find the image files
    image_files = dials_data("centroid_test_data").listdir("centroid*.cbf", sort=True)
    mask_filename = dials_data("centroid_test_data").join("mask.pickle")

    result = procrunner.run(
        [
            "dials.import",
            "mask=" + mask_filename.strpath,
            "output.experiments=experiments_with_mask.expt",
        ]
        + [f.strpath for f in image_files],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("experiments_with_mask.expt").check(file=1)

    from dxtbx.serialize import load

    experiments = load.experiment_list(
        tmpdir.join("experiments_with_mask.expt").strpath
    )
    assert (
        experiments[0].imageset.external_lookup.mask.filename == mask_filename.strpath
    )


def test_override_geometry(dials_data, tmpdir):
    # Find the image files
    image_files = dials_data("centroid_test_data").listdir("centroid*.cbf", sort=True)

    # Write a geometry phil file
    tmpdir.join("geometry.phil").write(
        """
      geometry {
        beam {
          wavelength = 2
          direction = (-1,0,0)
        }
        detector {
          panel {
            name = "New panel"
            type = "New type"
            pixel_size = 10,20
            image_size = 30,40
            trusted_range = 50,60
            thickness = 70
            material = "Si"
            fast_axis = -1,0,0
            slow_axis = 0,-1,0
            origin = 100,100,100
          }
        }
        goniometer {
          axes = 0,0,-1
          fixed_rotation = 0,1,2,3,4,5,6,7,8
          setting_rotation = 8,7,6,5,4,3,2,1,0
        }
        scan {
          image_range = 1,4
          oscillation = 1,2
        }
      }
  """
    )

    result = procrunner.run(
        ["dials.import", "geometry.phil", "output.experiments=override_geometry.expt"]
        + [f.strpath for f in image_files],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("override_geometry.expt").check(file=1)

    from dxtbx.serialize import load

    experiments = load.experiment_list(tmpdir.join("override_geometry.expt").strpath)
    imgset = experiments[0].imageset

    beam = imgset.get_beam()
    detector = imgset.get_detector()
    goniometer = imgset.get_goniometer()
    scan = imgset.get_scan()

    assert beam.get_wavelength() == 2
    assert beam.get_direction() == (-1, 0, 0)
    assert detector[0].get_name() == "New panel"
    assert detector[0].get_type() == "New type"
    assert detector[0].get_pixel_size() == (10, 20)
    assert detector[0].get_image_size() == (30, 40)
    assert detector[0].get_trusted_range() == (50, 60)
    assert detector[0].get_thickness() == 70
    assert detector[0].get_material() == "Si"
    assert detector[0].get_fast_axis() == (-1, 0, 0)
    assert detector[0].get_slow_axis() == (0, -1, 0)
    assert detector[0].get_origin() == (100, 100, 100)
    assert goniometer.get_rotation_axis_datum() == (0, 0, -1)
    assert goniometer.get_fixed_rotation() == (0, 1, 2, 3, 4, 5, 6, 7, 8)
    assert goniometer.get_setting_rotation() == (8, 7, 6, 5, 4, 3, 2, 1, 0)
    assert scan.get_image_range() == (1, 4)
    assert scan.get_oscillation() == (1, 2)


def tst_import_beam_centre(dials_data, run_in_tmpdir):
    # Find the image files
    image_files = dials_data("centroid_test_data").listdir("centroid*.cbf", sort=True)

    # provide mosflm beam centre to dials.import
    result = procrunner.run(
        [
            "dials.import",
            "mosflm_beam_centre=100,200",
            "output.experiments=mosflm_beam_centre.expt",
        ]
        + [f.strpath for f in image_files],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("mosflm_beam_centre.expt").check(file=1)

    from dxtbx.serialize import load

    experiments = load.experiment_list(tmpdir.join("mosflm_beam_centre.expt").strpath)
    imgset = experiments[0].imageset
    beam_centre = imgset.get_detector()[0].get_beam_centre(imgset.get_beam().get_s0())
    assert beam_centre == pytest.approx((200, 100))

    # provide an alternative models.expt to get geometry from
    result = procrunner.run(
        [
            "dials.import",
            "reference_geometry=mosflm_beam_centre.expt",
            "output.experiments=mosflm_beam_centre2.expt",
        ]
        + [f.strpath for f in image_files],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("mosflm_beam_centre2.expt").check(file=1)
    experiments = load.experiment_list(tmpdir.join("mosflm_beam_centre2.expt").strpath)
    imgset = experiments[0].imageset
    beam_centre = imgset.get_detector()[0].get_beam_centre(imgset.get_beam().get_s0())
    assert beam_centre == pytest.approx((200, 100))


def test_slow_fast_beam_centre(dials_regression, run_in_tmpdir):
    # test slow_fast_beam_centre with a multi-panel CS-PAD image
    impath = os.path.join(
        dials_regression,
        "image_examples",
        "LCLS_cspad_nexus",
        "idx-20130301060858401.cbf",
    )
    result = procrunner.run(
        [
            "dials.import",
            "slow_fast_beam_centre=134,42,18",
            "output.experiments=slow_fast_beam_centre.expt",
            impath,
        ]
    )
    assert not result.returncode and not result.stderr
    assert os.path.exists("slow_fast_beam_centre.expt")

    from dxtbx.serialize import load

    experiments = load.experiment_list("slow_fast_beam_centre.expt")
    imgset = experiments[0].imageset
    # beam centre on 18th panel
    s0 = imgset.get_beam().get_s0()
    beam_centre = imgset.get_detector()[18].get_beam_centre_px(s0)
    assert beam_centre == pytest.approx((42, 134))

    # check relative panel positions have not changed
    from scitbx import matrix

    o = matrix.col(imgset.get_detector()[0].get_origin())
    offsets = []
    for p in imgset.get_detector():
        intra_pnl = o - matrix.col(p.get_origin())
        offsets.append(intra_pnl.length())

    result = procrunner.run(
        ["dials.import", "output.experiments=reference.expt", impath]
    )
    assert not result.returncode and not result.stderr
    assert os.path.exists("reference.expt")

    ref_exp = load.experiment_list("reference.expt")
    ref_imset = ref_exp[0].imageset
    o = matrix.col(ref_imset.get_detector()[0].get_origin())
    ref_offsets = []
    for p in ref_imset.get_detector():
        intra_pnl = o - matrix.col(p.get_origin())
        ref_offsets.append(intra_pnl.length())
    assert offsets == pytest.approx(ref_offsets)


def test_from_image_files(dials_data, tmpdir):
    # Find the image files
    image_files = dials_data("centroid_test_data").listdir("centroid*.cbf", sort=True)

    # Import from the image files
    result = procrunner.run(
        ["dials.import", "output.experiments=imported.expt"]
        + [f.strpath for f in image_files],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode
    assert tmpdir.join("imported.expt").check(file=1)


def test_from_template(dials_data, tmpdir):
    # Find the image files
    template = dials_data("centroid_test_data").join("centroid_####.cbf")

    # Import from the image files
    result = procrunner.run(
        [
            "dials.import",
            "template=" + template.strpath,
            "output.experiments=imported.expt",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode
    assert tmpdir.join("imported.expt").check(file=1)


def test_extrapolate_scan(dials_data, tmpdir):
    # First image file
    image = dials_data("centroid_test_data").join("centroid_0001.cbf")

    result = procrunner.run(
        [
            "dials.import",
            image.strpath,
            "output.experiments=import_extrapolate.expt",
            "geometry.scan.image_range=1,900",
            "geometry.scan.extrapolate_scan=True",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode
    assert tmpdir.join("import_extrapolate.expt").check(file=1)
