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
def input_experiment_list(request, dials_data):
    filename = (dials_data(request.param) / "imported_experiments.json").strpath
    return load.experiment_list(filename)


def test_generate_mask(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.generate_mask",
            dials_data("centroid_test_data").join("experiments.json").strpath,
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result["exitcode"] and not result["stderr"]
    assert tmpdir.join("mask.pickle").check()


def test_generate_mask_with_untrusted_rectangle(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.generate_mask",
            dials_data("centroid_test_data").join("experiments.json").strpath,
            "output.mask=mask2.pickle",
            "output.experiments=masked_experiments.json",
            "untrusted.rectangle=100,200,100,200",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result["exitcode"] and not result["stderr"]
    assert tmpdir.join("mask2.pickle").check()
    assert tmpdir.join("masked_experiments.json").check()

    experiments = load.experiment_list(tmpdir.join("masked_experiments.json").strpath)
    imageset = experiments.imagesets()[0]
    assert imageset.external_lookup.mask.filename == tmpdir.join("mask2.pickle").strpath


def test_generate_mask_with_untrusted_circle(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.generate_mask",
            dials_data("centroid_test_data").join("experiments.json").strpath,
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
            dials_data("centroid_test_data").join("experiments.json").strpath,
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
            dials_data("centroid_test_data").join("experiments.json").strpath,
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
            dials_data("centroid_test_data").join("experiments.json").strpath,
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
            dials_data("centroid_test_data").join("experiments.json").strpath,
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


def test_generate_mask_function_with_untrusted_rectangle(input_experiment_list, tmpdir):
    params = phil_scope.extract()
    params.output.mask = tmpdir.join("mask4.pickle").strpath
    params.output.experiments = tmpdir.join("masked_experiments.json").strpath
    params.untrusted.rectangle = [100, 200, 100, 200]
    generate_mask(input_experiment_list, params)

    assert tmpdir.join("mask4.pickle").check() or all(
        [tmpdir.join("mask4_{:d}.pickle".format(i + 1)).check() for i in range(4)]
    )
    assert tmpdir.join("masked_experiments.json").check()

    experiments = load.experiment_list(tmpdir.join("masked_experiments.json").strpath)
    imageset = experiments.imagesets()[0]
    associated_masks = [
        tmpdir.join(f).strpath for f in ("mask4.pickle", "mask4_1.pickle")
    ]
    assert imageset.external_lookup.mask.filename in associated_masks


def test_apply_mask(dials_data, tmpdir):
    """
    Test applying a mask with the apply_mask parameter

    :param dials_data:  Use the centroid_test_data set from dials_data.
    :param tmpdir:  Perform the test in a temporary directory.
    """
    params = phil_scope.extract()
    input_expts = dials_data("centroid_test_data").join("experiments.json")
    input_mask = dials_data("centroid_test_data").join("mask.pickle")
    output_expts = tmpdir.join("masked_experiments.json")
    output_mask = tmpdir.join("mask5.pickle")
    params.apply_mask = [input_mask.strpath]
    params.output.mask = output_mask.strpath
    params.output.experiments = output_expts.strpath

    generate_mask(load.experiment_list(input_expts.strpath), params)

    # Test for the existence of the expected output files
    assert output_mask.check()
    assert output_expts.check()

    # Test that the mask hasn't been modified in any way
    with input_mask.open("rb") as old_maskfile:
        old_mask = pickle.load(old_maskfile)
    with output_mask.open("rb") as new_maskfile:
        new_mask = pickle.load(new_maskfile)
    assert old_mask == new_mask

    # Test that the mask has been associated with the imageset correctly
    imageset = load.experiment_list(output_expts.strpath).imagesets()[0]
    assert imageset.external_lookup.mask.filename == output_mask.strpath


def test_combine_masks(dials_data, tmpdir):
    """
    Test the use of apply_mask and other mask parameters to augment a mask

    :param dials_data:  Use the centroid_test_data set from dials_data.
    :param tmpdir:  Perform the test in a temporary directory.
    """
    params = phil_scope.extract()
    input_expts = dials_data("centroid_test_data").join("experiments.json")
    input_mask = dials_data("centroid_test_data").join("mask.pickle")
    output_expts = tmpdir.join("masked_experiments.json")
    output_mask = tmpdir.join("mask6.pickle")
    params.border = 100
    params.output.mask = output_mask.strpath
    params.output.experiments = output_expts.strpath

    # First create a simple 100px border mask
    generate_mask(load.experiment_list(input_expts.strpath), params)

    with output_mask.open("rb") as mask6_file:
        mask6 = pickle.load(mask6_file)

    # Then make another, this time augmenting an existing mask file
    output_mask = tmpdir.join("mask7.pickle")
    params.output.mask = output_mask.strpath
    params.apply_mask = [input_mask.strpath]
    generate_mask(load.experiment_list(input_expts.strpath), params)

    with output_mask.open("rb") as mask7_file:
        mask7 = pickle.load(mask7_file)
    with input_mask.open("rb") as mask_file:
        mask = pickle.load(mask_file)

    # Test that the augmentation of the existing mask with the border mask has simply
    # performed a boolean AND on each panel mask in mask.pickle with its counterpart
    # in mask6.pickle.
    combined_mask = tuple(old & new for old, new in zip(mask, mask6))
    assert mask7 == combined_mask
