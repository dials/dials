from __future__ import annotations

import random
from collections import namedtuple

import pytest

from dxtbx.serialize import load
from scitbx import matrix

from dials.algorithms.profile_model.gaussian_rs.transform import beam_vector_map

SequenceAndModel = namedtuple(
    "SequenceAndModel", ("sequence", "beam", "detector", "gonio", "scan")
)


@pytest.fixture
def sequence_and_model(dials_data):
    sequence = load.imageset(
        dials_data("centroid_test_data", pathlib=True) / "sweep.json"
    )
    return SequenceAndModel(
        sequence=sequence,
        beam=sequence.get_beam(),
        detector=sequence.get_detector(),
        gonio=sequence.get_goniometer(),
        scan=sequence.get_scan(),
    )


def test_at_corners(sequence_and_model):
    assert len(sequence_and_model.detector) == 1

    # The detector beam vectors
    ds1 = beam_vector_map(sequence_and_model.detector[0], sequence_and_model.beam, True)
    expected_size = sequence_and_model.detector[0].get_image_size()[::-1]
    expected_size = tuple(e + 1 for e in expected_size)
    assert ds1.all() == expected_size

    s0_length = matrix.col(sequence_and_model.beam.get_s0()).length()

    # Ensure a few random points are correct
    eps = 1e-7
    for k in range(1000):
        j = random.randint(0, ds1.all()[0] - 1)
        i = random.randint(0, ds1.all()[1] - 1)
        y = float(j)
        x = float(i)
        xyz = sequence_and_model.detector[0].get_pixel_lab_coord((x, y))
        s11 = matrix.col(xyz).normalize() * s0_length
        s12 = matrix.col(ds1[j, i])
        assert (s11 - s12).length() <= eps


def test_sub_division_at_corners(sequence_and_model):
    n_div = 2

    # The detector beam vectors
    ds1 = beam_vector_map(
        sequence_and_model.detector[0],
        sequence_and_model.beam,
        n_div,
        True,
    )
    expected_size = sequence_and_model.detector[0].get_image_size()[::-1]
    expected_size = tuple(e * n_div + 1 for e in expected_size)
    assert ds1.all() == expected_size

    s0_length = matrix.col(sequence_and_model.beam.get_s0()).length()

    # Ensure a few random points are correct
    eps = 1e-7
    for k in range(1000):
        j = random.randint(0, ds1.all()[0] - 1)
        i = random.randint(0, ds1.all()[1] - 1)
        y = float(j) / n_div
        x = float(i) / n_div
        xyz = sequence_and_model.detector[0].get_pixel_lab_coord((x, y))
        s11 = matrix.col(xyz).normalize() * s0_length
        s12 = matrix.col(ds1[j, i])
        assert (s11 - s12).length() <= eps


def test_sub_division_at_centres(sequence_and_model):
    n_div = 2

    # The detector beam vectors
    ds1 = beam_vector_map(
        sequence_and_model.detector[0],
        sequence_and_model.beam,
        n_div,
        False,
    )
    expected_size = sequence_and_model.detector[0].get_image_size()[::-1]
    expected_size = tuple(e * n_div for e in expected_size)
    assert ds1.all() == expected_size

    s0_length = matrix.col(sequence_and_model.beam.get_s0()).length()

    # Ensure a few random points are correct
    eps = 1e-7
    for k in range(1000):
        j = random.randint(0, ds1.all()[0] - 1)
        i = random.randint(0, ds1.all()[1] - 1)
        y = float(j + 0.5) / n_div
        x = float(i + 0.5) / n_div
        xyz = sequence_and_model.detector[0].get_pixel_lab_coord((x, y))
        s11 = matrix.col(xyz).normalize() * s0_length
        s12 = matrix.col(ds1[j, i])
        assert (s11 - s12).length() <= eps
