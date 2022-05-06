from __future__ import annotations

import pickle

import pytest

from dxtbx.model import ExperimentList
from libtbx import phil

from dials.command_line.dials_import import do_import
from dials.command_line.dials_import import phil_scope as import_phil_scope
from dials.command_line.generate_mask import generate_mask, phil_scope


@pytest.fixture(
    params=[
        {
            "directory": "centroid_test_data",
            "filename": "imported_experiments.json",
            "masks": ["pixels.mask"],
        },
        {
            "directory": "l_cysteine_dials_output",
            "filename": "imported.expt",
            "masks": ["pixels_%d.mask" % (i + 1) for i in range(4)],
        },
    ],
    ids=["One sequence", "Four sequences"],
)
def experiments_masks(request, dials_data):
    filename = (
        dials_data(request.param["directory"], pathlib=True) / request.param["filename"]
    )
    return ExperimentList.from_file(filename), request.param["masks"]


def test_generate_mask(experiments_masks, run_in_tmp_path):
    experiments, masks = experiments_masks

    params = phil_scope.fetch().extract()
    generate_mask(experiments, params)

    assert all(run_in_tmp_path.joinpath(mask).is_file() for mask in masks)


def test_generate_mask_with_untrusted_rectangle(experiments_masks, tmpdir):
    experiments, masks = experiments_masks

    params = phil_scope.fetch(
        phil.parse("untrusted.rectangle=100,200,100,200")
    ).extract()
    params.output.experiments = "masked.expt"
    with tmpdir.as_cwd():
        generate_mask(experiments, params)

    assert all(tmpdir.join(mask).check() for mask in masks)
    assert tmpdir.join("masked.expt").check()

    experiments = ExperimentList.from_file(tmpdir.join("masked.expt").strpath)
    imageset = experiments.imagesets()[0]
    assert imageset.external_lookup.mask.filename == tmpdir.join(masks[0]).strpath


def test_generate_mask_with_untrusted_circle(experiments_masks, tmpdir):
    experiments, masks = experiments_masks

    params = phil_scope.fetch(phil.parse("untrusted.circle=100,100,10")).extract()
    with tmpdir.as_cwd():
        generate_mask(experiments, params)

    assert all(tmpdir.join(mask).check() for mask in masks)


def test_generate_mask_with_resolution_range(experiments_masks, tmpdir):
    experiments, masks = experiments_masks

    params = phil_scope.fetch().extract()
    params.resolution_range = [(2, 3)]
    with tmpdir.as_cwd():
        generate_mask(experiments, params)

    assert all(tmpdir.join(mask).check() for mask in masks)


def test_generate_mask_with_d_min_d_max(experiments_masks, tmpdir):
    experiments, masks = experiments_masks

    params = phil_scope.fetch().extract()
    params.d_min = 3
    params.d_max = 2
    with tmpdir.as_cwd():
        generate_mask(experiments, params)

    assert all(tmpdir.join(mask).check() for mask in masks)


def test_generate_mask_with_ice_rings(experiments_masks, tmpdir):
    experiments, masks = experiments_masks

    params = phil_scope.fetch().extract()
    params.ice_rings.filter = True
    params.ice_rings.d_min = 2
    with tmpdir.as_cwd():
        generate_mask(experiments, params)

    assert all(tmpdir.join(mask).check() for mask in masks)


def test_generate_mask_with_untrusted_polygon_and_pixels(experiments_masks, tmpdir):
    experiments, masks = experiments_masks

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
        generate_mask(experiments, params)

    assert all(tmpdir.join(mask).check() for mask in masks)
    with tmpdir.join(masks[0]).open("rb") as fh:
        mask = pickle.load(fh)
    assert not mask[0][0, 0]
    assert not mask[0][1, 1]
    assert mask[0][0, 1]


def test_generate_mask_function_with_untrusted_rectangle(experiments_masks, tmpdir):
    experiments, masks = experiments_masks
    masks = [tmpdir.join(mask.replace("pixels", "pixels4")) for mask in masks]

    params = phil_scope.fetch().extract()
    params.output.mask = tmpdir.join("pixels4.mask").strpath
    params.output.experiments = tmpdir.join("masked.expt").strpath
    params.untrusted.rectangle = [100, 200, 100, 200]
    generate_mask(experiments, params)

    assert all(mask.check() for mask in masks)
    assert tmpdir.join("masked.expt").check()

    experiments = ExperimentList.from_file(tmpdir / "masked.expt")
    associated_masks = [
        imageset.external_lookup.mask.filename for imageset in experiments.imagesets()
    ]
    assert all(assoc_mask == mask for assoc_mask, mask in zip(associated_masks, masks))


def test_generate_mask_trusted_range(dials_data, tmpdir):
    # https://github.com/dials/dials/issues/978

    image_files = sorted(
        str(f) for f in dials_data("x4wide", pathlib=True).glob("*.cbf")
    )
    with tmpdir.as_cwd():
        # Import as usual
        do_import(
            ["output.experiments=no-overloads.expt"] + image_files,
            phil=import_phil_scope,
        )

        experiments = ExperimentList.from_file(tmpdir / "no-overloads.expt")
        params = phil_scope.fetch(
            phil.parse("untrusted.rectangle=100,200,100,200")
        ).extract()
        params.output.mask = "pixels1.mask"
        generate_mask(experiments, params)

        # Import with narrow trusted range to produce overloads
        do_import(
            ["trusted_range=-1,100", "output.experiments=overloads.expt"] + image_files,
            phil=import_phil_scope,
        )

        experiments = ExperimentList.from_file(tmpdir / "overloads.expt")
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


def test_generate_whole_panel_mask(experiments_masks, tmpdir):
    experiments, masks = experiments_masks

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
        generate_mask(experiments, params)

    assert all(tmpdir.join(mask).check() for mask in masks)
    with tmpdir.join(masks[0]).open("rb") as fh:
        mask = pickle.load(fh)
    assert mask[0].count(False) == len(mask[0])
