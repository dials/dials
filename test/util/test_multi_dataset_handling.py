"""
Tests for dials.util.multi_dataset_handling functions
"""
from __future__ import absolute_import, division, print_function

import pytest
from dials.array_family import flex
from dials.util.multi_dataset_handling import (
    assign_unique_identifiers,
    parse_multiple_datasets,
    select_datasets_on_ids,
)
from dxtbx.model import Experiment, ExperimentList


@pytest.fixture
def experiments():
    """Make a list of three empty experiments"""
    experiments = ExperimentList()
    experiments.append(Experiment())
    experiments.append(Experiment())
    experiments.append(Experiment())
    return experiments


@pytest.fixture
def experiments_024(experiments):
    experiments[0].identifier = "0"
    experiments[1].identifier = "2"
    experiments[2].identifier = "4"
    return experiments


@pytest.fixture
def reflections():
    """Make a list of three reflection tables"""
    rt1 = flex.reflection_table()
    rt1["id"] = flex.int([0, 0, 0])
    rt2 = flex.reflection_table()
    rt2["id"] = flex.int([1, 1])
    rt3 = flex.reflection_table()
    rt3["id"] = flex.int([4, 4])
    reflections = [rt1, rt2, rt3]
    return reflections


@pytest.fixture
def reflections_024(reflections):
    reflections[0].experiment_identifiers()[0] = "0"
    reflections[1].experiment_identifiers()[1] = "2"
    reflections[2].experiment_identifiers()[4] = "4"
    return reflections


def test_select_specific_datasets_using_id(experiments_024, reflections_024):
    use_datasets = ["0", "2"]
    experiments, refl = select_datasets_on_ids(
        experiments_024, reflections_024, use_datasets=use_datasets
    )
    assert len(experiments) == 2
    assert len(refl) == 2
    assert list(experiments.identifiers()) == ["0", "2"]


def test_exclude_specific_datasets_using_id(experiments_024, reflections_024):
    experiments, refl = select_datasets_on_ids(
        experiments_024, reflections_024, exclude_datasets=["0"]
    )
    assert len(refl) == 2
    assert list(experiments.identifiers()) == ["2", "4"]
    assert len(experiments) == 2


def test_raise_exception_when_selecting_and_excluding_datasets_at_same_time(
    experiments_024, reflections_024
):
    with pytest.raises(ValueError):
        experiments, refl = select_datasets_on_ids(
            experiments_024,
            reflections_024,
            use_datasets=["2", "4"],
            exclude_datasets=["0"],
        )


def test_raise_exception_when_excluding_non_existing_dataset(
    experiments_024, reflections_024
):
    with pytest.raises(ValueError):
        experiments, refl = select_datasets_on_ids(
            experiments_024, reflections_024, exclude_datasets=["1"]
        )


def test_selecting_everything_is_identity_function(experiments_024, reflections_024):
    exp, refl = select_datasets_on_ids(experiments_024, reflections_024)
    assert exp is experiments_024
    assert refl is reflections_024


def test_raise_exception_when_not_all_identifiers_set(experiments, reflections_024):
    experiments[0].identifier = "0"
    experiments[1].identifier = "2"
    with pytest.raises(ValueError):
        exp, refl = select_datasets_on_ids(
            experiments, reflections_024, use_datasets=["2"]
        )


def test_raise_exception_when_selecting_non_existing_dataset(
    experiments_024, reflections_024
):
    with pytest.raises(ValueError):
        exp, refl = select_datasets_on_ids(
            experiments_024, reflections_024, use_datasets=["3"]
        )


def test_correct_handling_with_multi_dataset_table(experiments_024):
    reflections = flex.reflection_table()
    reflections["id"] = flex.int([0, 1, 2])
    reflections.experiment_identifiers()[0] = "0"
    reflections.experiment_identifiers()[1] = "2"
    reflections.experiment_identifiers()[2] = "4"
    exp, refl = select_datasets_on_ids(
        experiments_024, [reflections], exclude_datasets=["2"]
    )
    assert list(refl[0].experiment_identifiers().values()) == ["0", "4"]
    assert list(refl[0]["id"]) == [0, 2]


def test_assignment_of_unique_identifiers_when_refl_table_ids_are_present(
    experiments, reflections
):
    assert list(experiments.identifiers()) == ["", "", ""]
    exp, rts = assign_unique_identifiers(experiments, reflections)
    expected_identifiers = ["0", "1", "2"]
    expected_ids = [0, 1, 2]
    # Check that identifiers are set in experiments and reflection table.
    assert list(exp.identifiers()) == expected_identifiers
    for i, refl in enumerate(rts):
        assert refl.experiment_identifiers()[i] == expected_identifiers[i]
        assert set(refl["id"]) == {i}
        assert i == expected_ids[i]


def test_assign_identifiers_where_none_are_set_but_refl_table_ids_have_duplicates(
    experiments, reflections
):
    reflections[2]["id"] = flex.int([0, 0])
    exp, rts = assign_unique_identifiers(experiments, reflections)
    expected_identifiers = ["0", "1", "2"]
    # Check that identifiers are set in experiments and reflection table.
    assert list(exp.identifiers()) == expected_identifiers
    expected_ids = [0, 1, 2]
    for i, refl in enumerate(rts):
        assert refl.experiment_identifiers()[i] == expected_identifiers[i]
        assert set(refl["id"]) == {i}
        assert i == expected_ids[i]


def test_raise_exception_when_existing_identifiers_are_inconsistent(
    experiments_024, reflections
):
    reflections[1].experiment_identifiers()[0] = "5"
    # should raise an assertion error for inconsistent identifiers
    with pytest.raises(ValueError):
        exp, rts = assign_unique_identifiers(experiments_024, reflections)


def test_cases_where_all_set_whether_reflection_table_is_split_or_not(
    experiments_024, reflections_024
):
    # should pass experiments back if identifiers all already set
    exp, rts = assign_unique_identifiers(experiments_024, reflections_024)
    expected_identifiers = ["0", "2", "4"]
    # Check that identifiers are set in experiments and reflection table.
    assert exp is experiments_024
    assert list(exp.identifiers()) == expected_identifiers
    expected_ids = [0, 1, 4]
    for i, refl in enumerate(rts):
        id_ = refl["id"][0]
        assert refl.experiment_identifiers()[id_] == expected_identifiers[i]
        assert set(refl["id"]) == {id_}
        assert id_ == expected_ids[i]


def test_raise_exception_if_unequal_experiments_and_reflections(experiments_024):
    reflections_multi = [flex.reflection_table()]
    reflections_multi[0]["id"] = flex.int([0, 1, 4])
    reflections_multi[0].experiment_identifiers()[0] = "0"
    reflections_multi[0].experiment_identifiers()[1] = "2"
    reflections_multi[0].experiment_identifiers()[4] = "4"
    del reflections_multi[0].experiment_identifiers()[4]
    del reflections_multi[0]["id"][2]
    with pytest.raises(ValueError):
        exp, rts = assign_unique_identifiers(experiments_024, reflections_multi)


def test_assigned_identifiers_are_kept_when_assigning_rest(experiments, reflections):
    # Now test that if some are set, these are maintained and unique ids are
    # set for the rest
    experiments[0].identifier = "1"
    reflections[0].experiment_identifiers()[0] = "1"
    exp, rts = assign_unique_identifiers(experiments, reflections)
    expected_identifiers = ["1", "0", "2"]
    assert list(exp.identifiers()) == expected_identifiers
    expected_ids = [0, 1, 2]
    for i, refl in enumerate(rts):
        id_ = refl["id"][0]
        assert refl.experiment_identifiers()[id_] == expected_identifiers[i]
        assert set(refl["id"]) == {id_}
        assert id_ == expected_ids[i]


def test_assigning_specified_identifiers(experiments, reflections):
    reflections[0].experiment_identifiers()[0] = "5"
    reflections[1].experiment_identifiers()[1] = "6"
    reflections[1].experiment_identifiers()[4] = "7"
    input_identifiers = ["0", "1", "10"]
    exp, rts = assign_unique_identifiers(experiments, reflections, input_identifiers)
    assert list(exp.identifiers()) == input_identifiers
    assert rts[0].experiment_identifiers()[0] == "0"
    assert rts[0]["id"][0] == 0
    assert rts[1].experiment_identifiers()[1] == "1"
    assert rts[1]["id"][0] == 1
    assert rts[2].experiment_identifiers()[2] == "10"
    assert rts[2]["id"][0] == 2

    # Test raises ValueError when wrong number of identifiers given
    with pytest.raises(ValueError):
        exp, rts = assign_unique_identifiers(experiments, reflections, ["0", "1"])


def test_parse_multiple_datasets():
    """Test the namesake function. This expects a list of reflection tables, and
    selects on the column named 'id'."""
    rt1 = flex.reflection_table()
    rt1["id"] = flex.int([0, 0, 0])
    rt1.experiment_identifiers()[0] = "0"
    rt2 = flex.reflection_table()
    rt2["id"] = flex.int([2, 2, 4, 4])
    rt2.experiment_identifiers()[2] = "2"
    rt2.experiment_identifiers()[4] = "4"
    single_tables = parse_multiple_datasets([rt2])
    assert len(single_tables) == 2
    assert list(set(single_tables[0]["id"])) == [2]
    assert list(set(single_tables[1]["id"])) == [4]
    single_tables = parse_multiple_datasets([rt1, rt2])
    assert list(set(single_tables[0]["id"])) == [0]
    assert list(set(single_tables[1]["id"])) == [2]
    assert list(set(single_tables[2]["id"])) == [4]
    assert len(single_tables) == 3
    single_tables = parse_multiple_datasets([rt1])
    assert len(single_tables) == 1
    assert list(set(single_tables[0]["id"])) == [0]

    # if a duplicate id is given, then this should be detected and new ids
    # determined for all datasets.
    rt3 = flex.reflection_table()
    rt3["id"] = flex.int([2, 2])
    rt3.experiment_identifiers()[2] = "5"
    single_tables = parse_multiple_datasets([rt1, rt2, rt3])
    assert len(single_tables) == 4
    assert list(set(single_tables[0]["id"])) == [0]
    assert single_tables[0].experiment_identifiers()[0] == "0"
    assert list(set(single_tables[1]["id"])) == [1]
    assert single_tables[1].experiment_identifiers()[1] == "2"
    assert list(set(single_tables[2]["id"])) == [2]
    assert single_tables[2].experiment_identifiers()[2] == "4"
    assert list(set(single_tables[3]["id"])) == [3]
    assert single_tables[3].experiment_identifiers()[3] == "5"
