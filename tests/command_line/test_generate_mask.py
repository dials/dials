from __future__ import annotations

import pickle

import pytest

from dxtbx.model import ExperimentList
from libtbx import phil

from dials.command_line.dials_import import do_import
from dials.command_line.dials_import import phil_scope as import_phil_scope
from dials.command_line.generate_mask import generate_mask, phil_scope, run


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
    filename = dials_data(request.param["directory"]) / request.param["filename"]
    return ExperimentList.from_file(filename), request.param["masks"]


def test_generate_mask(experiments_masks, run_in_tmp_path):
    experiments, masks = experiments_masks

    params = phil_scope.fetch().extract()
    generate_mask(experiments, params)

    assert all(run_in_tmp_path.joinpath(mask).is_file() for mask in masks)


def test_generate_mask_with_untrusted_rectangle(
    experiments_masks, tmp_path, monkeypatch
):
    experiments, masks = experiments_masks

    params = phil_scope.fetch(
        phil.parse("untrusted.rectangle=100,200,100,200")
    ).extract()
    params.output.experiments = "masked.expt"

    monkeypatch.chdir(tmp_path)
    generate_mask(experiments, params)

    assert all((tmp_path / mask).is_file() for mask in masks)
    assert (tmp_path / "masked.expt").is_file()

    experiments = ExperimentList.from_file(tmp_path / "masked.expt")
    imageset = experiments.imagesets()[0]
    assert imageset.external_lookup.mask.filename == str(tmp_path / masks[0])


def test_generate_mask_with_untrusted_circle(experiments_masks, tmp_path, monkeypatch):
    experiments, masks = experiments_masks

    params = phil_scope.fetch(phil.parse("untrusted.circle=100,100,10")).extract()
    monkeypatch.chdir(tmp_path)
    generate_mask(experiments, params)

    assert all((tmp_path / mask).is_file() for mask in masks)


def test_generate_mask_with_resolution_range(experiments_masks, tmp_path, monkeypatch):
    experiments, masks = experiments_masks

    params = phil_scope.fetch().extract()
    params.resolution_range = [(2, 3)]
    monkeypatch.chdir(tmp_path)
    generate_mask(experiments, params)

    assert all((tmp_path / mask).is_file() for mask in masks)


def test_generate_mask_with_d_min_d_max(experiments_masks, tmp_path, monkeypatch):
    experiments, masks = experiments_masks

    params = phil_scope.fetch().extract()
    params.d_min = 2
    params.d_max = 3
    monkeypatch.chdir(tmp_path)
    generate_mask(experiments, params)

    assert all((tmp_path / mask).is_file() for mask in masks)


def test_generate_mask_with_d_max_and_beam_at_pixel_centre(
    experiments_masks, tmp_path, monkeypatch
):
    # https://github.com/dials/dials/issues/2322

    experiments, masks = experiments_masks

    params = phil_scope.fetch().extract()
    params.d_max = 20

    # Modify experiment to put beam in the centre of a pixel
    beam = experiments[0].beam
    panel = experiments[0].detector[0]
    px_size = 0.1
    panel.set_pixel_size((px_size, px_size))  # ensure this is exact
    beam.set_s0((0, 0, -1))
    new_origin = (-1235.5 * px_size, 1279.5 * px_size, -190)
    panel.set_frame((1, 0, 0), (0, -1, 0), new_origin)
    assert (panel.get_beam_centre_px(beam.get_s0())) == (1235.5, 1279.5)

    params = phil_scope.fetch().extract()
    params.d_max = 10
    monkeypatch.chdir(tmp_path)
    generate_mask(experiments, params)

    assert all((tmp_path / mask).is_file() for mask in masks)


def test_generate_mask_with_ice_rings(experiments_masks, tmp_path, monkeypatch):
    experiments, masks = experiments_masks

    params = phil_scope.fetch().extract()
    params.ice_rings.filter = True
    params.ice_rings.d_min = 2
    monkeypatch.chdir(tmp_path)
    generate_mask(experiments, params)

    assert all((tmp_path / mask).is_file() for mask in masks)


def test_generate_mask_with_untrusted_polygon_and_pixels(
    experiments_masks, tmp_path, monkeypatch
):
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

    monkeypatch.chdir(tmp_path)
    generate_mask(experiments, params)

    assert all((tmp_path / mask).is_file() for mask in masks)
    with (tmp_path / masks[0]).open("rb") as fh:
        mask = pickle.load(fh)
    assert not mask[0][0, 0]
    assert not mask[0][1, 1]
    assert mask[0][0, 1]


def test_generate_mask_function_with_untrusted_rectangle(experiments_masks, tmp_path):
    experiments, masks = experiments_masks
    masks = [tmp_path / mask.replace("pixels", "pixels4") for mask in masks]

    params = phil_scope.fetch().extract()
    params.output.mask = str(tmp_path / "pixels4.mask")
    params.output.experiments = str(tmp_path / "masked.expt")
    params.untrusted.rectangle = [100, 200, 100, 200]
    generate_mask(experiments, params)

    assert all(mask.is_file() for mask in masks)
    assert (tmp_path / "masked.expt").is_file()

    experiments = ExperimentList.from_file(tmp_path / "masked.expt")
    associated_masks = [
        imageset.external_lookup.mask.filename for imageset in experiments.imagesets()
    ]
    assert all(
        assoc_mask == str(mask) for assoc_mask, mask in zip(associated_masks, masks)
    )


def test_generate_mask_trusted_range(dials_data, tmp_path, monkeypatch):
    # https://github.com/dials/dials/issues/978

    image_files = sorted(str(f) for f in dials_data("x4wide").glob("*.cbf"))
    monkeypatch.chdir(tmp_path)
    # Import as usual
    do_import(
        ["output.experiments=no-overloads.expt"] + image_files,
        phil=import_phil_scope,
    )

    experiments = ExperimentList.from_file(tmp_path / "no-overloads.expt")
    params = phil_scope.fetch(
        phil.parse("untrusted.rectangle=100,200,100,200")
    ).extract()
    params.output.mask = "pixels1.mask"
    generate_mask(experiments, params)

    # Import with narrow trusted range to produce overloads
    do_import(
        ["trusted_range=0,100", "output.experiments=overloads.expt"] + image_files,
        phil=import_phil_scope,
    )

    experiments = ExperimentList.from_file(tmp_path / "overloads.expt")
    params = phil_scope.fetch(
        phil.parse("untrusted.rectangle=100,200,100,200")
    ).extract()
    params.output.mask = "pixels2.mask"
    generate_mask(experiments, params)

    with (tmp_path / "pixels1.mask").open("rb") as fh:
        mask1 = pickle.load(fh)
    with (tmp_path / "pixels2.mask").open("rb") as fh:
        mask2 = pickle.load(fh)

    # Overloads should not be included in the mask
    assert (mask1[0] == mask2[0]).all_eq(True)


def test_generate_whole_panel_mask(experiments_masks, tmp_path, monkeypatch):
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

    monkeypatch.chdir(tmp_path)
    generate_mask(experiments, params)

    assert all((tmp_path / mask).is_file() for mask in masks)
    with (tmp_path / masks[0]).open("rb") as fh:
        mask = pickle.load(fh)
    assert mask[0].count(False) == len(mask[0])


def test_combine_masks(dials_data, run_in_tmp_path):
    path = dials_data("centroid_test_data")
    experiments_path = path / "imported_experiments.json"
    mask_path = path / "mask.pickle"
    experiments = ExperimentList.from_file(experiments_path)
    with (mask_path).open("rb") as fh:
        masks = [pickle.loads(fh.read(), encoding="bytes")]
    params = phil_scope.fetch().extract()

    # Combine with existing mask
    generate_mask(experiments, params, existing_masks=masks)
    assert run_in_tmp_path.joinpath("pixels.mask").is_file()

    # Combine only existing masks
    run(args=[str(mask_path), str(mask_path), "output.mask=pixels2.mask"])
    assert run_in_tmp_path.joinpath("pixels2.mask").is_file()
