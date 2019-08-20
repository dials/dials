from __future__ import absolute_import, division, print_function

import six.moves.cPickle as pickle
import os

import procrunner
import pytest

from dials.array_family import flex  # noqa: F401, import dependency


def test_find_spots_from_images(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.find_spots",
            "output.reflections=spotfinder.refl",
            "output.shoeboxes=True",
            "algorithm=dispersion",
        ]
        + [
            f.strpath for f in dials_data("centroid_test_data").listdir("centroid*.cbf")
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("spotfinder.refl").check(file=1)

    with tmpdir.join("spotfinder.refl").open("rb") as f:
        reflections = pickle.load(f)
    assert len(reflections) in range(653, 655)
    refl = reflections[0]
    assert refl["intensity.sum.value"] == pytest.approx(42)
    assert refl["bbox"] == pytest.approx((1398, 1400, 513, 515, 0, 1))
    assert refl["xyzobs.px.value"] == pytest.approx(
        (1399.1190476190477, 514.2142857142857, 0.5)
    )
    assert "shoebox" in reflections


def test_find_spots_with_resolution_filter(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.find_spots",
            "output.reflections=spotfinder.refl",
            "output.shoeboxes=False",
            "algorithm=dispersion",
            "filter.d_min=2",
            "filter.d_max=15",
        ]
        + [
            f.strpath for f in dials_data("centroid_test_data").listdir("centroid*.cbf")
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("spotfinder.refl").check(file=1)

    with tmpdir.join("spotfinder.refl").open("rb") as f:
        reflections = pickle.load(f)
    assert len(reflections) in range(467, 469)
    assert "shoebox" not in reflections


def test_find_spots_with_hot_mask(dials_data, tmpdir):
    # now write a hot mask
    result = procrunner.run(
        [
            "dials.find_spots",
            "write_hot_mask=True",
            "output.reflections=spotfinder.refl",
            "algorithm=dispersion",
            "output.shoeboxes=False",
        ]
        + [
            f.strpath for f in dials_data("centroid_test_data").listdir("centroid*.cbf")
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("spotfinder.refl").check(file=1)
    assert tmpdir.join("hot_mask_0.pickle").check(file=1)

    with tmpdir.join("spotfinder.refl").open("rb") as f:
        reflections = pickle.load(f)
    assert len(reflections) in range(653, 655)
    assert "shoebox" not in reflections

    with tmpdir.join("hot_mask_0.pickle").open("rb") as f:
        mask = pickle.load(f)
    assert len(mask) == 1
    assert mask[0].count(False) == 12


def test_find_spots_with_hot_mask_with_prefix(dials_data, tmpdir):
    # now write a hot mask
    result = procrunner.run(
        [
            "dials.find_spots",
            "write_hot_mask=True",
            "hot_mask_prefix=my_hot_mask",
            "output.reflections=spotfinder.refl",
            "output.shoeboxes=False",
            "algorithm=dispersion",
        ]
        + [
            f.strpath for f in dials_data("centroid_test_data").listdir("centroid*.cbf")
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("spotfinder.refl").check(file=1)
    assert tmpdir.join("my_hot_mask_0.pickle").check(file=1)

    with tmpdir.join("spotfinder.refl").open("rb") as f:
        reflections = pickle.load(f)
    assert len(reflections) in range(653, 655)
    assert "shoebox" not in reflections
    with tmpdir.join("my_hot_mask_0.pickle").open("rb") as f:
        mask = pickle.load(f)
    assert len(mask) == 1
    assert mask[0].count(False) == 12


def test_find_spots_with_generous_parameters(dials_data, tmpdir):
    # now with more generous parameters
    result = procrunner.run(
        [
            "dials.find_spots",
            "min_spot_size=3",
            "max_separation=3",
            "output.reflections=spotfinder.refl",
            "algorithm=dispersion",
        ]
        + [
            f.strpath for f in dials_data("centroid_test_data").listdir("centroid*.cbf")
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("spotfinder.refl").check(file=1)

    with tmpdir.join("spotfinder.refl").open("rb") as f:
        reflections = pickle.load(f)
    assert len(reflections) in range(678, 680)


def test_find_spots_with_user_defined_mask(dials_data, tmpdir):
    # Now with a user defined mask
    result = procrunner.run(
        [
            "dials.find_spots",
            "output.reflections=spotfinder.refl",
            "output.shoeboxes=True",
            "algorithm=dispersion",
            "lookup.mask="
            + dials_data("centroid_test_data").join("mask.pickle").strpath,
        ]
        + [
            f.strpath for f in dials_data("centroid_test_data").listdir("centroid*.cbf")
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("spotfinder.refl").check(file=1)

    with tmpdir.join("spotfinder.refl").open("rb") as f:
        reflections = pickle.load(f)

    from dxtbx.model.experiment_list import ExperimentListFactory

    experiments = ExperimentListFactory.from_json_file(
        dials_data("centroid_test_data").join("experiments.json").strpath
    )
    assert len(experiments) == 1
    imageset = experiments.imagesets()[0]
    detector = imageset.get_detector()
    beam = imageset.get_beam()
    for x, y, z in reflections["xyzobs.px.value"]:
        d = detector[0].get_resolution_at_pixel(beam.get_s0(), (x, y))
        assert d >= 3


def test_find_spots_with_user_defined_region(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.find_spots",
            "output.reflections=spotfinder.refl",
            "output.shoeboxes=True",
            "region_of_interest=800,1200,800,1200",
        ]
        + [
            f.strpath for f in dials_data("centroid_test_data").listdir("centroid*.cbf")
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("spotfinder.refl").check(file=1)

    with tmpdir.join("spotfinder.refl").open("rb") as f:
        reflections = pickle.load(f)
    x, y, z = reflections["xyzobs.px.value"].parts()
    assert x.all_ge(800)
    assert y.all_ge(800)
    assert x.all_lt(1200)
    assert y.all_lt(1200)


def test_find_spots_with_xfel_stills(dials_regression, tmpdir):
    # now with XFEL stills
    result = procrunner.run(
        [
            "dials.find_spots",
            os.path.join(
                dials_regression,
                "spotfinding_test_data",
                "idx-s00-20131106040302615.cbf",
            ),
            "output.reflections=spotfinder.refl",
            "algorithm=dispersion",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("spotfinder.refl").check(file=1)

    with tmpdir.join("spotfinder.refl").open("rb") as f:
        reflections = pickle.load(f)
    assert len(reflections) == 2643
