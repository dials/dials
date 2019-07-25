from __future__ import absolute_import, division, print_function

import procrunner


def test_spots_xds(tmpdir):
    xds_input = "SPOT.XDS"
    output_pickle = "spot.refl"

    tmpdir.join(xds_input).write(
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
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join(output_pickle).check(file=1)

    from dials.array_family import flex

    reflections = flex.reflection_table.from_file(tmpdir.join(output_pickle).strpath)
    assert len(reflections) == 5

    tmpdir.join(xds_input).remove()
    assert not tmpdir.join(xds_input).check()

    # now test we can export it again
    result = procrunner.run(
        ["dials.export", "format=xds", output_pickle], working_directory=tmpdir.strpath
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("xds", "SPOT.XDS").check(file=1)

    txt = tmpdir.join("xds", "SPOT.XDS").read()
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


def test_export_xds(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.find_spots",
            dials_data("centroid_test_data").join("experiments.json").strpath,
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("strong.refl").check(file=1)

    result = procrunner.run(
        [
            "dials.export",
            "format=xds",
            dials_data("centroid_test_data").join("experiments.json").strpath,
            "strong.refl",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("xds", "XDS.INP").check(file=1)
    assert tmpdir.join("xds", "XPARM.XDS").check(file=1)
    assert tmpdir.join("xds", "SPOT.XDS").check(file=1)

    tmpdir.join("xds", "XDS.INP").remove()
    tmpdir.join("xds", "XPARM.XDS").remove()
    assert not tmpdir.join("xds", "XDS.INP").check(file=1)
    assert not tmpdir.join("xds", "XPARM.XDS").check(file=1)
    result = procrunner.run(
        [
            "dials.export",
            "format=xds",
            dials_data("centroid_test_data").join("experiments.json").strpath,
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("xds", "XDS.INP").check(file=1)
    assert tmpdir.join("xds", "XPARM.XDS").check(file=1)
