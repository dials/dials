"""
Tests for the functions in dials.util.options
"""

from __future__ import annotations

from unittest.mock import Mock

import pytest

from dxtbx.model import Experiment, ExperimentList

from dials.util import Sorry
from dials.util.options import (
    ArgumentParser,
    flatten_reflections,
    reflections_and_experiments_from_files,
)

from . import mock_reflection_file_object, mock_two_reflection_file_object


def test_cannot_read_headerless_h5(dials_data):
    data_h5 = dials_data("vmxi_thaumatin", pathlib=True) / "image_15799_data_000001.h5"
    parser = ArgumentParser(read_experiments_from_images=True)
    with pytest.raises(Sorry):
        parser.parse_args([str(data_h5)])


def test_flatten_experiments_updating_id_values():
    """Test the correct handling of duplicate table id values.

    Note that this function does not have the ability to update the
    experiment string identifier, only ensure that the table id values
    do not clash (it is not possible even to load multiple experiments
    with the same identifier).
    """
    # Test the case of two single reflection tables.
    file_list = [mock_reflection_file_object(id_=0), mock_reflection_file_object(id_=0)]
    rs = flatten_reflections(file_list)
    assert rs[0] is file_list[0].data
    assert list(rs[0]["id"]) == [-1, 0, 0]
    assert list(rs[0].experiment_identifiers().keys()) == [0]
    assert list(rs[0].experiment_identifiers().values()) == ["0"]
    assert rs[1] is file_list[1].data
    assert list(rs[1]["id"]) == [-1, 1, 1]
    assert list(rs[1].experiment_identifiers().keys()) == [1]
    assert list(rs[1].experiment_identifiers().values()) == ["0"]

    # Now test the case where one reflection table contains two experiments
    file_list = [mock_two_reflection_file_object(), mock_reflection_file_object(id_=0)]
    rs = flatten_reflections(file_list)
    assert rs[0] is file_list[0].data
    assert list(rs[0]["id"]) == [-1, 0, 0, 1, 1]
    assert list(rs[0].experiment_identifiers().keys()) == [0, 1]
    assert list(rs[0].experiment_identifiers().values()) == ["0", "2"]
    assert rs[1] is file_list[1].data
    assert list(rs[1]["id"]) == [-1, 2, 2]
    assert list(rs[1].experiment_identifiers().keys()) == [2]
    assert list(rs[1].experiment_identifiers().values()) == ["0"]

    file_list = [
        mock_reflection_file_object(id_=0),
        mock_two_reflection_file_object(ids=[1, 2]),
    ]
    rs = flatten_reflections(file_list)
    assert rs[0] is file_list[0].data
    assert list(rs[0]["id"]) == [-1, 0, 0]
    assert list(rs[0].experiment_identifiers().keys()) == [0]
    assert list(rs[0].experiment_identifiers().values()) == ["0"]
    assert rs[1] is file_list[1].data
    assert list(rs[1]["id"]) == [-1, 1, 1, 2, 2]
    assert list(rs[1].experiment_identifiers().keys()) == [1, 2]
    assert list(rs[1].experiment_identifiers().values()) == ["1", "2"]


def test_reflections_and_experiments_from_files():
    """Test correct extracting of reflections and experiments."""
    # Test when input reflections order matches the experiments order
    refl_file_list = [
        mock_two_reflection_file_object(ids=[0, 1]),
        mock_reflection_file_object(id_=2),
    ]

    def mock_exp_obj(id_=0):
        """Make a mock experiments file object."""
        exp = Mock()
        exp.data = ExperimentList()
        exp.data.append(Experiment(identifier=str(id_)))
        return exp

    exp_file_list = [mock_exp_obj(id_=i) for i in [0, 1, 2]]

    refls, expts = reflections_and_experiments_from_files(refl_file_list, exp_file_list)
    assert refls[0] is refl_file_list[0].data
    assert refls[1] is refl_file_list[1].data
    assert expts[0].identifier == "0"
    assert expts[1].identifier == "1"
    assert expts[2].identifier == "2"

    # Test when input reflections order does not match experiments order.
    refl_file_list = [
        mock_reflection_file_object(id_=2),
        mock_two_reflection_file_object(ids=[0, 1]),
    ]
    refls, expts = reflections_and_experiments_from_files(refl_file_list, exp_file_list)
    assert refls[0] is refl_file_list[1].data
    assert refls[1] is refl_file_list[0].data
    assert expts[0].identifier == "0"
    assert expts[1].identifier == "1"
    assert expts[2].identifier == "2"
