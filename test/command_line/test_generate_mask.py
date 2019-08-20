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
    assert not result.returncode and not result.stderr
    assert tmpdir.join("pixels.mask").check()


def test_generate_mask_with_untrusted_rectangle(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.generate_mask",
            dials_data("centroid_test_data").join("experiments.json").strpath,
            "output.mask=pixels2.mask",
            "output.experiments=masked.expt",
            "untrusted.rectangle=100,200,100,200",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("pixels2.mask").check()
    assert tmpdir.join("masked.expt").check()

    experiments = load.experiment_list(tmpdir.join("masked.expt").strpath)
    imageset = experiments.imagesets()[0]
    assert imageset.external_lookup.mask.filename == tmpdir.join("pixels2.mask").strpath


def test_generate_mask_with_untrusted_circle(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.generate_mask",
            dials_data("centroid_test_data").join("experiments.json").strpath,
            "output.mask=pixels3.mask",
            "untrusted.circle=100,100,10",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("pixels3.mask").check()


def test_generate_mask_with_resolution_range(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.generate_mask",
            dials_data("centroid_test_data").join("experiments.json").strpath,
            "output.mask=pixels4.mask",
            "resolution_range=2,3",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("pixels4.mask").check()


def test_generate_mask_with_d_min_d_max(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.generate_mask",
            dials_data("centroid_test_data").join("experiments.json").strpath,
            "output.mask=pixels5.mask",
            "d_min=3",
            "d_max=2",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("pixels5.mask").check()


def test_generate_mask_with_ice_rings(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.generate_mask",
            dials_data("centroid_test_data").join("experiments.json").strpath,
            "output.mask=pixels6.mask",
            "ice_rings{filter=True;d_min=2}",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("pixels6.mask").check()


def test_generate_mask_with_untrusted_polygon_and_pixels(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.generate_mask",
            dials_data("centroid_test_data").join("experiments.json").strpath,
            "output.mask=pixels3.mask",
            "untrusted.polygon=100,100,100,200,200,200,200,100",
            "untrusted.pixel=0,0",
            "untrusted.pixel=1,1",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("pixels3.mask").check()
    with tmpdir.join("pixels3.mask").open("rb") as fh:
        mask = pickle.load(fh)
    assert not mask[0][0, 0]
    assert not mask[0][1, 1]
    assert mask[0][0, 1]


def test_generate_mask_function_with_untrusted_rectangle(input_experiment_list, tmpdir):
    params = phil_scope.extract()
    params.output.mask = tmpdir.join("pixels4.mask").strpath
    params.output.experiments = tmpdir.join("masked.expt").strpath
    params.untrusted.rectangle = [100, 200, 100, 200]
    generate_mask(input_experiment_list, params)

    assert tmpdir.join("pixels4.mask").check() or all(
        [tmpdir.join("pixels4_{:d}.mask".format(i + 1)).check() for i in range(4)]
    )
    assert tmpdir.join("masked.expt").check()

    experiments = load.experiment_list(tmpdir.join("masked.expt").strpath)
    imageset = experiments.imagesets()[0]
    associated_masks = [
        tmpdir.join(f).strpath for f in ("pixels4.mask", "pixels4_1.mask")
    ]
    assert imageset.external_lookup.mask.filename in associated_masks
