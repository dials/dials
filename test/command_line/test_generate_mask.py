from __future__ import absolute_import, division, print_function

import six.moves.cPickle as pickle
import pytest

from dials.command_line.generate_mask import generate_mask, phil_scope
from dials.command_line.dials_import import Script as ImportScript
from dials.command_line.dials_import import phil_scope as import_phil_scope
from dxtbx.serialize import load
from libtbx import phil


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


def test_generate_mask(input_experiment_list, tmpdir):
    params = phil_scope.fetch().extract()
    with tmpdir.as_cwd():
        generate_mask(input_experiment_list, params)

    assert tmpdir.join("pixels.mask").check()


def test_generate_mask_with_untrusted_rectangle(input_experiment_list, tmpdir):
    params = phil_scope.fetch(
        phil.parse("untrusted.rectangle=100,200,100,200")
    ).extract()
    params.output.experiments = "masked.expt"
    with tmpdir.as_cwd():
        generate_mask(input_experiment_list, params)

    assert tmpdir.join("pixels.mask").check()
    assert tmpdir.join("masked.expt").check()

    experiments = load.experiment_list(tmpdir.join("masked.expt").strpath)
    imageset = experiments.imagesets()[0]
    assert imageset.external_lookup.mask.filename == tmpdir.join("pixels.mask").strpath


def test_generate_mask_with_untrusted_circle(input_experiment_list, tmpdir):
    params = phil_scope.fetch(phil.parse("untrusted.circle=100,100,10")).extract()
    with tmpdir.as_cwd():
        generate_mask(input_experiment_list, params)

    assert tmpdir.join("pixels.mask").check()


def test_generate_mask_with_resolution_range(input_experiment_list, tmpdir):
    params = phil_scope.fetch().extract()
    params.resolution_range = [(2, 3)]
    with tmpdir.as_cwd():
        generate_mask(input_experiment_list, params)

    assert tmpdir.join("pixels.mask").check()


def test_generate_mask_with_d_min_d_max(input_experiment_list, tmpdir):
    params = phil_scope.fetch().extract()
    params.d_min = 3
    params.d_max = 2
    with tmpdir.as_cwd():
        generate_mask(input_experiment_list, params)

    assert tmpdir.join("pixels.mask").check()


def test_generate_mask_with_ice_rings(input_experiment_list, tmpdir):
    params = phil_scope.fetch().extract()
    params.ice_rings.filter = True
    params.ice_rings.d_min = 2
    with tmpdir.as_cwd():
        generate_mask(input_experiment_list, params)

    assert tmpdir.join("pixels.mask").check()


def test_generate_mask_with_untrusted_polygon_and_pixels(input_experiment_list, tmpdir):
    params = phil_scope.fetch(
        phil.parse(
            """
untrusted {
  polygon = 100 100 100 200 200 200 200 100
}
untrusted {
  pixel = 0 0
}
untrusted {
  pixel = 1 1
}"""
        )
    ).extract()

    with tmpdir.as_cwd():
        generate_mask(input_experiment_list, params)

    assert tmpdir.join("pixels.mask").check()
    with tmpdir.join("pixels.mask").open("rb") as fh:
        mask = pickle.load(fh)
    assert not mask[0][0, 0]
    assert not mask[0][1, 1]
    assert mask[0][0, 1]


def test_generate_mask_function_with_untrusted_rectangle(input_experiment_list, tmpdir):
    params = phil_scope.fetch().extract()
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

        experiments = load.experiment_list(tmpdir.join("no-overloads.expt").strpath)
        params = phil_scope.fetch(
            phil.parse("untrusted.rectangle=100,200,100,200")
        ).extract()
        params.output.mask = "pixels1.mask"
        generate_mask(experiments, params)

        # Import with narrow trusted range to produce overloads
        import_script = ImportScript(import_phil_scope)
        import_script.run(
            ["trusted_range=-1,100", "output.experiments=overloads.expt"] + image_files
        )

        experiments = load.experiment_list(tmpdir.join("overloads.expt").strpath)
        params = phil_scope.fetch(
            phil.parse("untrusted.rectangle=100,200,100,200")
        ).extract()
        params.output.mask = "pixels2.mask"
        generate_mask(experiments, params)

    with tmpdir.join("pixels1.mask").open("rb") as fh:
        mask1 = pickle.load(fh)
    with tmpdir.join("pixels2.mask").open("rb") as fh:
        mask2 = pickle.load(fh)

    # Overloads should not be included in the mask
    assert (mask1[0] == mask2[0]).all_eq(True)


def test_generate_whole_panel_mask(input_experiment_list, tmpdir):
    params = phil_scope.fetch(
        phil.parse(
            """
untrusted {
  panel = 0
}
"""
        )
    ).extract()

    with tmpdir.as_cwd():
        generate_mask(input_experiment_list, params)

    assert tmpdir.join("pixels.mask").check()
    with tmpdir.join("pixels.mask").open("rb") as fh:
        mask = pickle.load(fh)
    assert mask[0].count(False) == len(mask[0])
