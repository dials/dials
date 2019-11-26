from __future__ import absolute_import, division, print_function

import six.moves.cPickle as pickle
import pytest

from dials.command_line.generate_mask import generate_mask, run, phil_scope
from dials.command_line.dials_import import Script as ImportScript
from dials.command_line.dials_import import phil_scope as import_phil_scope
from dxtbx.serialize import load


@pytest.fixture(
    params=[
        "centroid_test_data",
        pytest.param(
            "l_cysteine_dials_output",
            marks=pytest.mark.skip(
                reason="test depends on https://github.com/dials/data/pull/40"
            ),
        ),
    ],
    ids=["One sequence", "Four sequences"],
)
def input_experiment_list(request, dials_data):
    filename = (dials_data(request.param) / "imported_experiments.json").strpath
    return load.experiment_list(filename)


def test_generate_mask(dials_data, tmpdir):
    with tmpdir.as_cwd():
        run(
            phil_scope,
            [dials_data("centroid_test_data").join("experiments.json").strpath],
        )

    assert tmpdir.join("pixels.mask").check()


def test_generate_mask_with_untrusted_rectangle(dials_data, tmpdir):
    with tmpdir.as_cwd():
        run(
            phil_scope,
            [
                dials_data("centroid_test_data").join("experiments.json").strpath,
                "output.mask=pixels2.mask",
                "output.experiments=masked.expt",
                "untrusted.rectangle=100,200,100,200",
            ],
        )

    assert tmpdir.join("pixels2.mask").check()
    assert tmpdir.join("masked.expt").check()

    experiments = load.experiment_list(tmpdir.join("masked.expt").strpath)
    imageset = experiments.imagesets()[0]
    assert imageset.external_lookup.mask.filename == tmpdir.join("pixels2.mask").strpath


def test_generate_mask_with_untrusted_circle(dials_data, tmpdir):
    with tmpdir.as_cwd():
        run(
            phil_scope,
            [
                dials_data("centroid_test_data").join("experiments.json").strpath,
                "output.mask=pixels3.mask",
                "untrusted.circle=100,100,10",
            ],
        )

    assert tmpdir.join("pixels3.mask").check()


def test_generate_mask_with_resolution_range(dials_data, tmpdir):
    with tmpdir.as_cwd():
        run(
            phil_scope,
            [
                dials_data("centroid_test_data").join("experiments.json").strpath,
                "output.mask=pixels4.mask",
                "resolution_range=2,3",
            ],
        )

    assert tmpdir.join("pixels4.mask").check()


def test_generate_mask_with_d_min_d_max(dials_data, tmpdir):
    with tmpdir.as_cwd():
        run(
            phil_scope,
            [
                dials_data("centroid_test_data").join("experiments.json").strpath,
                "output.mask=pixels5.mask",
                "d_min=3",
                "d_max=2",
            ],
        )

    assert tmpdir.join("pixels5.mask").check()


def test_generate_mask_with_ice_rings(dials_data, tmpdir):
    with tmpdir.as_cwd():
        run(
            phil_scope,
            [
                dials_data("centroid_test_data").join("experiments.json").strpath,
                "output.mask=pixels6.mask",
                "ice_rings{filter=True;d_min=2}",
            ],
        )

    assert tmpdir.join("pixels6.mask").check()


def test_generate_mask_with_untrusted_polygon_and_pixels(dials_data, tmpdir):
    with tmpdir.as_cwd():
        run(
            phil_scope,
            [
                dials_data("centroid_test_data").join("experiments.json").strpath,
                "output.mask=pixels3.mask",
                "untrusted.polygon=100,100,100,200,200,200,200,100",
                "untrusted.pixel=0,0",
                "untrusted.pixel=1,1",
            ],
        )

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
        tmpdir.join("pixels4_{:d}.mask".format(i + 1)).check() for i in range(4)
    )
    assert tmpdir.join("masked.expt").check()

    experiments = load.experiment_list(tmpdir.join("masked.expt").strpath)
    imageset = experiments.imagesets()[0]
    associated_masks = [
        tmpdir.join(f).strpath for f in ("pixels4.mask", "pixels4_1.mask")
    ]
    assert imageset.external_lookup.mask.filename in associated_masks


def test_generate_mask_trusted_range(dials_data, tmpdir):
    # https://github.com/dials/dials/issues/978

    image_files = [f.strpath for f in dials_data("x4wide").listdir("*.cbf", sort=True)]
    with tmpdir.as_cwd():
        # Import as usual
        import_script = ImportScript(import_phil_scope)
        import_script.run(["output.experiments=no-overloads.expt"] + image_files)

        run(
            phil_scope,
            [
                "no-overloads.expt",
                "output.mask=pixels1.mask",
                "untrusted.rectangle=100,200,100,200",
            ],
        )

        # Import with narrow trusted range to produce overloads
        import_script = ImportScript(import_phil_scope)
        import_script.run(
            ["trusted_range=-1,100", "output.experiments=overloads.expt"] + image_files
        )

        run(
            phil_scope,
            [
                "overloads.expt",
                "output.mask=pixels2.mask",
                "untrusted.rectangle=100,200,100,200",
            ],
        )

    with tmpdir.join("pixels1.mask").open("rb") as fh:
        mask1 = pickle.load(fh)
    with tmpdir.join("pixels2.mask").open("rb") as fh:
        mask2 = pickle.load(fh)

    # Overloads should not be included in the mask
    assert (mask1[0] == mask2[0]).all_eq(True)
