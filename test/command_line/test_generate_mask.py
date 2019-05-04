from __future__ import absolute_import, division, print_function

import six.moves.cPickle as pickle
import procrunner
import pytest

from dials.command_line.generate_mask import generate_mask, phil_scope
from dxtbx.serialize import load


@pytest.fixture(
    params=[
        "centroid_test_data",
        pytest.param("l_cysteine_dials_output", marks=pytest.mark.xfail),
    ],
    ids=["One sweep", "Four sweeps"],
)
def input_datablock(request, dials_data):
    filename = (dials_data(request.param) / "datablock.json").strpath
    return load.datablock(filename)[0]


def test_generate_mask(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.generate_mask",
            dials_data("centroid_test_data").join("datablock.json").strpath,
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result["exitcode"] and not result["stderr"]
    assert tmpdir.join("mask.pickle").check()


def test_generate_mask_with_untrusted_rectangle(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.generate_mask",
            dials_data("centroid_test_data").join("datablock.json").strpath,
            "output.mask=mask2.pickle",
            "output.datablock=masked_datablock.json",
            "untrusted.rectangle=100,200,100,200",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result["exitcode"] and not result["stderr"]
    assert tmpdir.join("mask2.pickle").check()
    assert tmpdir.join("masked_datablock.json").check()

    datablock = load.datablock(tmpdir.join("masked_datablock.json").strpath)[0]
    imageset = datablock.extract_imagesets()[0]
    assert imageset.external_lookup.mask.filename == tmpdir.join("mask2.pickle").strpath


def test_generate_mask_with_untrusted_circle(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.generate_mask",
            dials_data("centroid_test_data").join("datablock.json").strpath,
            "output.mask=mask3.pickle",
            "untrusted.circle=100,100,10",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result["exitcode"] and not result["stderr"]
    assert tmpdir.join("mask3.pickle").check()


def test_generate_mask_with_resolution_range(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.generate_mask",
            dials_data("centroid_test_data").join("datablock.json").strpath,
            "output.mask=mask4.pickle",
            "resolution_range=2,3",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result["exitcode"] and not result["stderr"]
    assert tmpdir.join("mask4.pickle").check()


def test_generate_mask_with_d_min_d_max(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.generate_mask",
            dials_data("centroid_test_data").join("datablock.json").strpath,
            "output.mask=mask5.pickle",
            "d_min=3",
            "d_max=2",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result["exitcode"] and not result["stderr"]
    assert tmpdir.join("mask5.pickle").check()


def test_generate_mask_with_ice_rings(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.generate_mask",
            dials_data("centroid_test_data").join("datablock.json").strpath,
            "output.mask=mask6.pickle",
            "ice_rings{filter=True;d_min=2}",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result["exitcode"] and not result["stderr"]
    assert tmpdir.join("mask6.pickle").check()


def test_generate_mask_with_untrusted_polygon_and_pixels(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.generate_mask",
            dials_data("centroid_test_data").join("datablock.json").strpath,
            "output.mask=mask3.pickle",
            "untrusted.polygon=100,100,100,200,200,200,200,100",
            "untrusted.pixel=0,0",
            "untrusted.pixel=1,1",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result["exitcode"] and not result["stderr"]
    assert tmpdir.join("mask3.pickle").check()
    with tmpdir.join("mask3.pickle").open("rb") as fh:
        mask = pickle.load(fh)
    assert not mask[0][0, 0]
    assert not mask[0][1, 1]
    assert mask[0][0, 1]


def test_generate_mask_function_with_untrusted_rectangle(input_datablock, tmpdir):
    params = phil_scope.extract()
    params.output.mask = tmpdir.join("mask4.pickle").strpath
    params.output.datablock = tmpdir.join("masked_datablock.json").strpath
    params.untrusted.rectangle = [100, 200, 100, 200]
    generate_mask(input_datablock, params)

    assert tmpdir.join("mask4.pickle").check() or all(
        [tmpdir.join("mask4_{:d}.pickle".format(i + 1)).check() for i in range(4)]
    )
    assert tmpdir.join("masked_datablock.json").check()

    datablock = load.experiment_list(tmpdir.join("masked_datablock.json").strpath)
    imageset = datablock.imagesets()[0]
    associated_masks = [
        tmpdir.join(f).strpath for f in ("mask4.pickle", "mask4_1.pickle")
    ]
    assert imageset.external_lookup.mask.filename in associated_masks
