from __future__ import annotations

import os

import procrunner


def test_export_single_bitmap(dials_data, tmp_path):
    result = procrunner.run(
        [
            "dials.export_bitmaps",
            dials_data("centroid_test_data", pathlib=True) / "centroid_0001.cbf",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("image0001.png").is_file()


def test_export_multiple_bitmaps(dials_data, tmp_path):
    result = procrunner.run(
        [
            "dials.export_bitmaps",
            dials_data("centroid_test_data", pathlib=True) / "experiments.json",
            "prefix=variance_",
            "binning=2",
            "display=variance",
            "colour_scheme=inverse_greyscale",
            "brightness=25",
            "kernel_size=5,5",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr

    for i in range(1, 8):
        assert tmp_path.joinpath(f"variance_000{i}.png").is_file()


def test_export_bitmap_with_prefix_and_no_padding(dials_data, tmp_path):
    result = procrunner.run(
        [
            "dials.export_bitmaps",
            dials_data("centroid_test_data", pathlib=True) / "centroid_0001.cbf",
            "prefix=img_",
            "padding=0",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("img_1.png").is_file()


def test_export_bitmap_with_prefix_and_extra_padding(dials_data, tmp_path):
    result = procrunner.run(
        [
            "dials.export_bitmaps",
            dials_data("centroid_test_data", pathlib=True) / "centroid_0001.cbf",
            "prefix=img_",
            "padding=5",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("img_00001.png").is_file()


def test_export_bitmap_with_specified_output_filename(dials_data, tmp_path):
    result = procrunner.run(
        [
            "dials.export_bitmaps",
            dials_data("centroid_test_data", pathlib=True) / "centroid_0001.cbf",
            "output.file=kittens.png",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("kittens.png").is_file()


def test_export_multiple_bitmaps_with_specified_output_filename_fails(
    dials_data, tmp_path
):
    # setting output filename not allowed with >1 image
    result = procrunner.run(
        [
            "dials.export_bitmaps",
            dials_data("centroid_test_data", pathlib=True) / "experiments.json",
            "output.file=kittens.png",
        ],
        working_directory=tmp_path,
    )
    assert result.returncode


def test_export_still_image(dials_regression, tmp_path):
    image = os.path.join(
        dials_regression, "image_examples", "DLS_I24_stills", "still_0001.cbf"
    )

    result = procrunner.run(["dials.export_bitmaps", image], working_directory=tmp_path)
    assert not result.returncode and not result.stderr

    assert tmp_path.joinpath("image0001.png").is_file()


def test_export_multi_panel(dials_regression, tmp_path):
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
            working_directory=tmp_path,
        )
        assert not result.returncode and not result.stderr
        assert tmp_path.joinpath(f"binning_{binning}_0001.png").is_file()


def test_export_restricted_multiimage(dials_data, tmp_path):
    "Test exporting a subset of an imageset"
    result = procrunner.run(
        [
            "dials.export_bitmaps",
            dials_data("centroid_test_data", pathlib=True) / "experiments.json",
            "imageset_index=2",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert [f.name for f in tmp_path.glob("*.png")] == [
        "image0002.png"
    ], "Only one image expected"
