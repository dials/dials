from __future__ import absolute_import, division, print_function

import procrunner


def test(dials_data, tmpdir):
    input_filename = dials_data("centroid_test_data").join("datablock.json").strpath

    result = procrunner.run(
        ["dials.estimate_gain", "input.experiments=" + input_filename],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert "Estimated gain: 1.0" in result["stdout"]
