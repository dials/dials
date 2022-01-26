"""Test for new experiment identifier features"""
from __future__ import annotations

from dxtbx.model import Experiment, ExperimentList

from dials.array_family import flex


def test_selection_identifier_propagation():
    """
    If experiment idenfitiers are set, then when selecting on a reflection
    table, copy across any for the same id.

    If no identifiers are set, then nothing should happen."""

    # Test when identifiers & id column set.
    refl = flex.reflection_table()
    refl["id"] = flex.int([0, 0, 0, 1, 1, 1, 2, 2, 2])
    refl.experiment_identifiers()[0] = "0"
    refl.experiment_identifiers()[1] = "1"
    refl.experiment_identifiers()[2] = "2"

    new_refl = refl.select(refl["id"] == 1)
    assert list(new_refl.experiment_identifiers().keys()) == [1]
    assert list(new_refl.experiment_identifiers().values()) == ["1"]

    # test with no identifiers set, but there is an id column
    refl = flex.reflection_table()
    refl["id"] = flex.int([0, 0, 0, 1, 1, 1, 2, 2, 2])
    new_refl = refl.select(refl["id"] == 0)
    assert new_refl.ncols() == 1
    assert list(new_refl.experiment_identifiers().keys()) == []
    assert list(new_refl.experiment_identifiers().values()) == []

    # Test selections still work when nothing id related is set.
    refl = flex.reflection_table()
    refl["x"] = flex.int([0, 0, 0, 1, 1, 1, 2, 2, 2])
    new_refl = refl.select(refl["x"] == 0)
    assert new_refl.ncols() == 1
    assert list(new_refl.experiment_identifiers().keys()) == []
    assert list(new_refl.experiment_identifiers().values()) == []


def test_select_using_experiment():
    """
    want to be able to do:
    for exp in experiments:
        refl = reflections.select(exp)

    rather than:
    for i, _ in enumerate(experiments):
        refl = reflections.select(refl["id"] == i)
    """
    exp1 = Experiment(identifier="1")
    exp2 = Experiment(identifier="2")

    refl = flex.reflection_table()
    refl["id"] = flex.int([0, 0, 0, 1, 1, 1])
    refl["idx"] = flex.int([0, 1, 2, 3, 4, 5])
    refl.experiment_identifiers()[0] = "0"
    refl.experiment_identifiers()[1] = "1"

    sel_refl = refl.select(exp1)
    assert list(sel_refl.experiment_identifiers().keys()) == [1]
    assert list(sel_refl.experiment_identifiers().values()) == ["1"]
    assert list(sel_refl["idx"]) == [3, 4, 5]

    # Test returns nothing if no match
    sel_refl = refl.select(exp2)
    assert list(sel_refl.experiment_identifiers().keys()) == []
    assert list(sel_refl.experiment_identifiers().values()) == []
    assert not sel_refl

    # Test returns nothing if no identifiers set. Thought is that this is
    # more useful than raising a dxtbx assert error, as it allows handling within
    # python in case identifiers not set.
    refl = flex.reflection_table()
    refl["idx"] = flex.int([0, 1, 2, 3, 4, 5])
    sel_refl = refl.select(exp1)
    assert list(sel_refl.experiment_identifiers().keys()) == []
    assert list(sel_refl.experiment_identifiers().values()) == []
    assert not sel_refl


def test_select_using_experiments():
    """Test selecting using an experimentlist

    e.g refl = reflections.select(experiments) to select a subset of the
    data from an experimentlist."""

    refl = flex.reflection_table()
    refl["id"] = flex.int([0, 0, 0, 1, 1, 1, 2, 2, 2])
    refl["idx"] = flex.int([0, 1, 2, 3, 4, 5, 6, 7, 8])
    refl.experiment_identifiers()[0] = "0"
    refl.experiment_identifiers()[1] = "1"
    refl.experiment_identifiers()[2] = "2"
    experiments = ExperimentList(
        [
            Experiment(identifier="0"),
            Experiment(identifier="2"),
            Experiment(identifier="3"),
        ]
    )

    sel_refl = refl.select(experiments)
    assert list(sel_refl.experiment_identifiers().keys()) == [0, 2]
    assert list(sel_refl.experiment_identifiers().values()) == ["0", "2"]
    assert list(sel_refl["idx"]) == [0, 1, 2, 6, 7, 8]
