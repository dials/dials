from __future__ import absolute_import, division, print_function

import random

import pytest

from scitbx import matrix

from dxtbx.serialize import load

from dials.algorithms.profile_model.gaussian_rs.transform import beam_vector_map


@pytest.fixture
def sequence_and_model(dials_data):
    class Test(object):
        pass

    storage_class = Test()

    storage_class.sequence = load.imageset(
        dials_data("centroid_test_data").join("sweep.json").strpath
    )

    # Get the models
    storage_class.beam = storage_class.sequence.get_beam()
    storage_class.detector = storage_class.sequence.get_detector()
    storage_class.gonio = storage_class.sequence.get_goniometer()
    storage_class.scan = storage_class.sequence.get_scan()

    return storage_class


def test_at_corners(sequence_and_model):
    assert len(sequence_and_model.detector) == 1

    # The detector beam vectors
    ds1 = beam_vector_map(sequence_and_model.detector[0], sequence_and_model.beam, True)
    expected_size = sequence_and_model.detector[0].get_image_size()[::-1]
    expected_size = tuple([e + 1 for e in expected_size])
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
    sequence_and_model.n_div = 2

    # The detector beam vectors
    ds1 = beam_vector_map(
        sequence_and_model.detector[0],
        sequence_and_model.beam,
        sequence_and_model.n_div,
        True,
    )
    expected_size = sequence_and_model.detector[0].get_image_size()[::-1]
    expected_size = tuple([e * sequence_and_model.n_div + 1 for e in expected_size])
    assert ds1.all() == expected_size

    s0_length = matrix.col(sequence_and_model.beam.get_s0()).length()

    # Ensure a few random points are correct
    eps = 1e-7
    for k in range(1000):
        j = random.randint(0, ds1.all()[0] - 1)
        i = random.randint(0, ds1.all()[1] - 1)
        y = float(j) / sequence_and_model.n_div
        x = float(i) / sequence_and_model.n_div
        xyz = sequence_and_model.detector[0].get_pixel_lab_coord((x, y))
        s11 = matrix.col(xyz).normalize() * s0_length
        s12 = matrix.col(ds1[j, i])
        assert (s11 - s12).length() <= eps


def test_sub_division_at_centres(sequence_and_model):
    sequence_and_model.n_div = 2

    # The detector beam vectors
    ds1 = beam_vector_map(
        sequence_and_model.detector[0],
        sequence_and_model.beam,
        sequence_and_model.n_div,
        False,
    )
    expected_size = sequence_and_model.detector[0].get_image_size()[::-1]
    expected_size = tuple([e * sequence_and_model.n_div for e in expected_size])
    assert ds1.all() == expected_size

    s0_length = matrix.col(sequence_and_model.beam.get_s0()).length()

    # Ensure a few random points are correct
    eps = 1e-7
    for k in range(1000):
        j = random.randint(0, ds1.all()[0] - 1)
        i = random.randint(0, ds1.all()[1] - 1)
        y = float(j + 0.5) / sequence_and_model.n_div
        x = float(i + 0.5) / sequence_and_model.n_div
        xyz = sequence_and_model.detector[0].get_pixel_lab_coord((x, y))
        s11 = matrix.col(xyz).normalize() * s0_length
        s12 = matrix.col(ds1[j, i])
        assert (s11 - s12).length() <= eps
