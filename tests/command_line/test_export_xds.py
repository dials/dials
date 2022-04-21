from __future__ import annotations

import procrunner

from dials.array_family import flex


def test_spots_xds(tmp_path):
    xds_input = "SPOT.XDS"
    output_pickle = "spot.refl"

    tmp_path.joinpath(xds_input).write_text(
        """\
 2411.40 1000.70 25.00 16384. 0 0 0
 1328.60 2170.40 20.57 7326. 0 0 0
 177.56 2191.30 24.94 6779. 0 0 0
 1231.34 1220.04 24.99 1952. 0 0 0
 1227.07 1230.56 24.81 349. 0 0 0
 1341.63 1243.25 5.64 321. 2 -2 11
 1125.23 1197.72 12.14 231. -1 2 -10
 1317.52 1171.59 19.28 120. 6 -4 6
 1260.25 1300.55 13.67 116. -4 2 6
 1090.27 1199.47 41.49 114. -2 3 -13
"""
    )

    result = procrunner.run(
        [
            "dials.import_xds",
            xds_input,  # xparm_file,
            "input.method=reflections",
            "output.filename=" + output_pickle,
            "remove_invalid=True",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath(output_pickle).is_file()

    reflections = flex.reflection_table.from_file(tmp_path / output_pickle)
    assert len(reflections) == 5

    tmp_path.joinpath(xds_input).unlink()
    assert not tmp_path.joinpath(xds_input).exists()

    # now test we can export it again
    result = procrunner.run(
        ["dials.export", "format=xds", output_pickle], working_directory=tmp_path
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("xds", "SPOT.XDS").is_file()

    txt = tmp_path.joinpath("xds", "SPOT.XDS").read_text()
    assert [line.strip() for line in txt.split("\n")] == [
        line.strip()
        for line in """\
 1341.63 1243.25 5.64 321.00  2 -2 11
 1125.23 1197.72 12.14 231.00  -1 2 -10
 1317.52 1171.59 19.28 120.00  6 -4 6
 1260.25 1300.55 13.67 116.00  -4 2 6
 1090.27 1199.47 41.49 114.00  -2 3 -13
""".split(
            "\n"
        )
    ]


def test_export_xds(dials_data, tmp_path):
    experiment = dials_data("centroid_test_data", pathlib=True) / "experiments.json"
    result = procrunner.run(
        ["dials.find_spots", "nproc=1", experiment], working_directory=tmp_path
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("strong.refl").is_file()

    result = procrunner.run(
        [
            "dials.export",
            "format=xds",
            experiment,
            "strong.refl",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("xds", "XDS.INP").is_file()
    assert tmp_path.joinpath("xds", "XPARM.XDS").is_file()
    assert tmp_path.joinpath("xds", "SPOT.XDS").is_file()

    tmp_path.joinpath("xds", "XDS.INP").unlink()
    tmp_path.joinpath("xds", "XPARM.XDS").unlink()
    assert not tmp_path.joinpath("xds", "XDS.INP").is_file()
    assert not tmp_path.joinpath("xds", "XPARM.XDS").is_file()
    result = procrunner.run(
        [
            "dials.export",
            "format=xds",
            experiment,
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("xds", "XDS.INP").is_file()
    assert tmp_path.joinpath("xds", "XPARM.XDS").is_file()


def test_export_imported_experiments(dials_data, tmp_path):
    result = procrunner.run(
        [
            "dials.export",
            "format=xds",
            dials_data("centroid_test_data", pathlib=True)
            / "imported_experiments.json",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("xds", "XDS.INP").is_file()
    assert not tmp_path.joinpath("xds", "XPARM.XDS").is_file()
    assert not tmp_path.joinpath("xds", "SPOT.XDS").is_file()
