from __future__ import annotations

import os
import pickle
import shutil

import procrunner
import pytest

from dxtbx.model.experiment_list import ExperimentListFactory

import dials.command_line.find_spots
from dials.array_family import flex


def _check_expected_results(reflections):
    """Check expected results for standard processing."""
    assert len(reflections) in range(653, 655)
    refl = reflections[0]
    assert refl["intensity.sum.value"] == pytest.approx(42)
    assert refl["bbox"] == pytest.approx((1398, 1400, 513, 515, 0, 1))
    assert refl["xyzobs.px.value"] == pytest.approx(
        (1399.1190476190477, 514.2142857142857, 0.5)
    )
    assert "shoebox" in reflections


def test_find_spots_from_images(dials_data, tmp_path):
    result = procrunner.run(
        [
            "dials.find_spots",
            "nproc=1",
            "output.reflections=spotfinder.refl",
            "output.shoeboxes=True",
            "algorithm=dispersion",
        ]
        + list(dials_data("centroid_test_data", pathlib=True).glob("centroid*.cbf")),
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "spotfinder.refl").is_file()

    reflections = flex.reflection_table.from_file(tmp_path / "spotfinder.refl")
    _check_expected_results(reflections)
    # No identifiers set if just running on images and not outputting experiments
    assert not reflections.experiment_identifiers().values(), list(
        reflections.experiment_identifiers().values()
    )


def test_find_spots_from_images_override_maximum(dials_data, tmp_path):
    result = procrunner.run(
        [
            "dials.find_spots",
            "nproc=1",
            "maximum_trusted_value=100",
            "output.reflections=spotfinder.refl",
            "output.shoeboxes=True",
            "algorithm=dispersion",
        ]
        + list(dials_data("centroid_test_data", pathlib=True).glob("centroid*.cbf")),
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "spotfinder.refl").is_file()

    reflections = flex.reflection_table.from_file(tmp_path / "spotfinder.refl")
    sbox = reflections["shoebox"]
    for s in sbox:
        assert flex.max(s.data) <= 100


def test_find_spots_from_zero_indexed_cbf(dials_data, tmp_path):
    one_indexed_cbf = (
        dials_data("centroid_test_data", pathlib=True) / "centroid_0001.cbf"
    )
    zero_indexed_cbf = tmp_path / "centroid_0000.cbf"
    shutil.copy(one_indexed_cbf, zero_indexed_cbf)

    result = procrunner.run(
        ["dials.find_spots", "nproc=1", zero_indexed_cbf], working_directory=tmp_path
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "strong.refl").is_file()
    assert b"Saved 0 reflections to" not in result.stdout, "No spots found on 0000.cbf"


def test_find_spots_from_images_output_experiments(dials_data, tmp_path):
    result = procrunner.run(
        [
            "dials.find_spots",
            "nproc=1",
            "output.reflections=spotfinder.refl",
            "output.shoeboxes=True",
            "algorithm=dispersion",
            "output.experiments=spotfinder.expt",
        ]
        + list(dials_data("centroid_test_data", pathlib=True).glob("centroid*.cbf")),
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "spotfinder.refl").is_file()

    reflections = flex.reflection_table.from_file(tmp_path / "spotfinder.refl")
    _check_expected_results(reflections)
    # Identifiers set if experiments are output
    assert reflections.experiment_identifiers().values(), list(
        reflections.experiment_identifiers().values()
    )


def test_find_spots_from_imported_experiments(dials_data, tmp_path):
    """First run import to generate an imported.expt and use this."""
    _ = procrunner.run(
        ["dials.import"]
        + list(dials_data("centroid_test_data", pathlib=True).glob("centroid*.cbf")),
        working_directory=tmp_path,
    )

    result = procrunner.run(
        [
            "dials.find_spots",
            "nproc=1",
            tmp_path / "imported.expt",
            "output.reflections=spotfinder.refl",
            "output.shoeboxes=True",
            "algorithm=dispersion",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "spotfinder.refl").is_file()

    reflections = flex.reflection_table.from_file(tmp_path / "spotfinder.refl")
    _check_expected_results(reflections)
    # Identifiers set if just running on images
    assert len(reflections.experiment_identifiers().values()) == 1, list(
        reflections.experiment_identifiers().values()
    )


def test_find_spots_from_imported_as_grid(dials_data, tmp_path):
    """First run import to generate an imported.expt and use this."""
    _ = procrunner.run(
        ["dials.import", "oscillation=0,0"]
        + list(dials_data("centroid_test_data", pathlib=True).glob("centroid*.cbf")),
        working_directory=tmp_path,
    )

    result = procrunner.run(
        [
            "dials.find_spots",
            "nproc=1",
            tmp_path / "imported.expt",
            "output.reflections=spotfinder.refl",
            "output.shoeboxes=True",
            "algorithm=dispersion",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "spotfinder.refl").is_file()

    reflections = flex.reflection_table.from_file(tmp_path / "spotfinder.refl")
    assert len(set(reflections["id"])) == 9, len(set(reflections["id"]))


def test_find_spots_with_resolution_filter(dials_data, tmp_path):
    result = procrunner.run(
        [
            "dials.find_spots",
            "nproc=1",
            "output.reflections=spotfinder.refl",
            "output.shoeboxes=False",
            "algorithm=dispersion",
            "filter.d_min=2",
            "filter.d_max=15",
        ]
        + list(dials_data("centroid_test_data", pathlib=True).glob("centroid*.cbf")),
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "spotfinder.refl").is_file()

    reflections = flex.reflection_table.from_file(tmp_path / "spotfinder.refl")
    assert len(reflections) in range(467, 469)
    assert "shoebox" not in reflections


def test_find_spots_with_hot_mask(dials_data, tmp_path):
    # now write a hot mask
    result = procrunner.run(
        [
            "dials.find_spots",
            "nproc=1",
            "write_hot_mask=True",
            "output.reflections=spotfinder.refl",
            "algorithm=dispersion",
            "output.shoeboxes=False",
        ]
        + list(dials_data("centroid_test_data", pathlib=True).glob("centroid*.cbf")),
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "spotfinder.refl").is_file()
    assert (tmp_path / "hot_mask_0.pickle").is_file()

    reflections = flex.reflection_table.from_file(tmp_path / "spotfinder.refl")
    assert len(reflections) in range(653, 655)
    assert "shoebox" not in reflections

    with (tmp_path / "hot_mask_0.pickle").open("rb") as f:
        mask = pickle.load(f)
    assert len(mask) == 1
    assert mask[0].count(False) == 12


def test_find_spots_with_hot_mask_with_prefix(dials_data, tmp_path):
    # now write a hot mask
    result = procrunner.run(
        [
            "dials.find_spots",
            "nproc=1",
            "write_hot_mask=True",
            "hot_mask_prefix=my_hot_mask",
            "output.reflections=spotfinder.refl",
            "output.shoeboxes=False",
            "algorithm=dispersion",
        ]
        + list(dials_data("centroid_test_data", pathlib=True).glob("centroid*.cbf")),
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "spotfinder.refl").is_file()
    assert (tmp_path / "my_hot_mask_0.pickle").is_file()

    reflections = flex.reflection_table.from_file(tmp_path / "spotfinder.refl")
    assert len(reflections) in range(653, 655)
    assert "shoebox" not in reflections
    with (tmp_path / "my_hot_mask_0.pickle").open("rb") as f:
        mask = pickle.load(f)
    assert len(mask) == 1
    assert mask[0].count(False) == 12


def test_find_spots_with_generous_parameters(dials_data, tmp_path):
    # now with more generous parameters
    result = procrunner.run(
        [
            "dials.find_spots",
            "nproc=1",
            "min_spot_size=3",
            "max_separation=3",
            "output.reflections=spotfinder.refl",
            "algorithm=dispersion",
        ]
        + list(dials_data("centroid_test_data", pathlib=True).glob("centroid*.cbf")),
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "spotfinder.refl").is_file()

    reflections = flex.reflection_table.from_file(tmp_path / "spotfinder.refl")
    assert len(reflections) in range(678, 680)


def test_find_spots_with_user_defined_mask(dials_data, tmp_path):
    # Now with a user defined mask
    mask_pickle = dials_data("centroid_test_data", pathlib=True) / "mask.pickle"
    result = procrunner.run(
        [
            "dials.find_spots",
            "nproc=1",
            "output.reflections=spotfinder.refl",
            "output.shoeboxes=True",
            "algorithm=dispersion",
            f"lookup.mask={mask_pickle}",
        ]
        + list(dials_data("centroid_test_data", pathlib=True).glob("centroid*.cbf")),
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "spotfinder.refl").is_file()

    reflections = flex.reflection_table.from_file(tmp_path / "spotfinder.refl")
    experiments = ExperimentListFactory.from_json_file(
        dials_data("centroid_test_data", pathlib=True) / "experiments.json"
    )
    assert len(experiments) == 1
    imageset = experiments.imagesets()[0]
    detector = imageset.get_detector()
    beam = imageset.get_beam()
    for x, y, z in reflections["xyzobs.px.value"]:
        d = detector[0].get_resolution_at_pixel(beam.get_s0(), (x, y))
        assert d >= 3


def test_find_spots_with_user_defined_region(dials_data, tmp_path):
    result = procrunner.run(
        [
            "dials.find_spots",
            "nproc=1",
            "output.reflections=spotfinder.refl",
            "output.shoeboxes=True",
            "region_of_interest=800,1200,800,1200",
        ]
        + list(dials_data("centroid_test_data", pathlib=True).glob("centroid*.cbf")),
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "spotfinder.refl").is_file()

    reflections = flex.reflection_table.from_file(tmp_path / "spotfinder.refl")
    x, y, z = reflections["xyzobs.px.value"].parts()
    assert x.all_ge(800)
    assert y.all_ge(800)
    assert x.all_lt(1200)
    assert y.all_lt(1200)


def test_find_spots_with_xfel_stills(dials_regression, tmp_path):
    # now with XFEL stills
    result = procrunner.run(
        [
            "dials.find_spots",
            "nproc=1",
            os.path.join(
                dials_regression,
                "spotfinding_test_data",
                "idx-s00-20131106040302615.cbf",
            ),
            "output.reflections=spotfinder.refl",
            "algorithm=dispersion",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "spotfinder.refl").is_file()

    reflections = flex.reflection_table.from_file(tmp_path / "spotfinder.refl")
    assert len(reflections) == 2643


def test_find_spots_with_per_image_statistics(dials_data, tmp_path):
    result = procrunner.run(
        ["dials.find_spots", "nproc=1", "per_image_statistics=True"]
        + list(dials_data("centroid_test_data", pathlib=True).glob("centroid*.cbf")),
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "strong.refl").is_file()
    assert b"Number of centroids per image for imageset 0:" in result.stdout
    assert (
        b"|   image |   #spots |   #spots_no_ice |   total_intensity |" in result.stdout
    )


@pytest.mark.parametrize(
    "blur,expected_nref", [("None", 559), ("narrow", 721), ("wide", 739)]
)
def test_find_spots_radial_profile(dials_data, blur, expected_nref, run_in_tmp_path):
    reflections = dials.command_line.find_spots.run(
        [
            "nproc=1",
            "threshold.algorithm=radial_profile",
            f"blur={blur}",
        ]
        + [
            os.fspath(f)
            for f in dials_data("centroid_test_data", pathlib=True).glob(
                "centroid*.cbf"
            )
        ],
        return_results=True,
    )
    assert len(reflections) == expected_nref
