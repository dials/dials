"""
Tests for the functions in dials.util.options
"""
from __future__ import absolute_import, division, print_function

from dials.util.options import flatten_reflections, flatten_experiments, OptionParser
from dials.test.util import mock_reflection_file_object, mock_two_reflection_file_object


def test_can_read_headerless_h5_and_no_detector_is_present(dials_data):
    data_h5 = dials_data("vmxi_thaumatin").join("image_15799_data_000001.h5").strpath
    parser = OptionParser(read_experiments=True, read_experiments_from_images=True)
    params, options = parser.parse_args([data_h5])
    experiments = flatten_experiments(params.input.experiments)
    assert len(experiments) == 1
    assert not experiments[0].detector


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
