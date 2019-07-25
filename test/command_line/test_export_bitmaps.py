from __future__ import absolute_import, division, print_function

import os
import procrunner


def test_export_single_bitmap(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.export_bitmaps",
            dials_data("centroid_test_data").join("centroid_0001.cbf").strpath,
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("image0001.png").check(file=1)


def test_export_multiple_bitmaps(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.export_bitmaps",
            dials_data("centroid_test_data").join("experiments.json").strpath,
            "prefix=variance_",
            "binning=2",
            "display=variance",
            "colour_scheme=inverse_greyscale",
            "brightness=25",
            "kernel_size=5,5",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr

    for i in range(1, 8):
        assert tmpdir.join("variance_000%i.png" % i).check(file=1)


def test_export_bitmap_with_prefix_and_no_padding(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.export_bitmaps",
            dials_data("centroid_test_data").join("centroid_0001.cbf").strpath,
            "prefix=img_",
            "padding=0",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("img_1.png").check(file=1)


def test_export_bitmap_with_prefix_and_extra_padding(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.export_bitmaps",
            dials_data("centroid_test_data").join("centroid_0001.cbf").strpath,
            "prefix=img_",
            "padding=5",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("img_00001.png").check(file=1)


def test_export_bitmap_with_specified_output_filename(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.export_bitmaps",
            dials_data("centroid_test_data").join("centroid_0001.cbf").strpath,
            "output_file=kittens.png",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("kittens.png").check(file=1)


def test_export_multiple_bitmaps_with_specified_output_filename_fails(
    dials_data, tmpdir
):
    # setting output filename not allowed with >1 image
    result = procrunner.run(
        [
            "dials.export_bitmaps",
            dials_data("centroid_test_data").join("experiments.json").strpath,
            "output_file=kittens.png",
        ],
        working_directory=tmpdir.strpath,
    )
    assert result.returncode


def test_export_still_image(dials_regression, tmpdir):
    image = os.path.join(
        dials_regression, "image_examples", "DLS_I24_stills", "still_0001.cbf"
    )

    result = procrunner.run(
        ["dials.export_bitmaps", image], working_directory=tmpdir.strpath
    )
    assert not result.returncode and not result.stderr

    assert tmpdir.join("image0001.png").check(file=1)


def test_export_multi_panel(dials_regression, tmpdir):
    image = os.path.join(
        dials_regression, "image_examples", "DLS_I23", "germ_13KeV_0001.cbf"
    )

    for binning in (1, 4):
        result = procrunner.run(
            [
                "dials.export_bitmaps",
                image,
                "binning=%i" % binning,
                "prefix=binning_%i_" % binning,
            ],
            working_directory=tmpdir.strpath,
        )
        assert not result.returncode and not result.stderr

        assert tmpdir.join("binning_%i_0001.png" % binning).check(file=1)


def test_export_restricted_multiimage(dials_data, tmpdir):
    "Test exporting a subset of an imageset"
    result = procrunner.run(
        [
            "dials.export_bitmaps",
            dials_data("centroid_test_data").join("experiments.json").strpath,
            "imageset_index=2",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert [f.basename for f in tmpdir.listdir("*.png")] == [
        "image0002.png"
    ], "Only one image expected"
