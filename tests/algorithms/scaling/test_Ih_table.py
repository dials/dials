from __future__ import annotations

from unittest.mock import Mock

import numpy as np
import pandas as pd
import pytest

from cctbx.sgtbx import space_group, uctbx
from dxtbx import flumpy
from scitbx import sparse

from dials.algorithms.scaling.Ih_table import IhTable, IhTableBlock, map_indices_to_asu
from dials.array_family import flex


@pytest.fixture()
def large_reflection_table():
    """Create a reflection table."""
    return generate_refl_1()


@pytest.fixture()
def small_reflection_table():
    """Create a small reflection table."""
    return generate_refl_2()


@pytest.fixture(scope="module")
def test_sg():
    """Create a space group object."""
    return space_group("C 2y")


def mock_error_model():
    """Mock error model."""
    em = Mock()
    em.update_variances.return_value = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
    return em


def refl_table_to_df(table):
    df = pd.DataFrame()
    for col in table.keys():
        if col != "asu_miller_index" and col != "miller_index":
            df[col] = flumpy.to_numpy(table[col])
    return df


def generate_refl_1():
    """Generate a test reflection table. Note that the variance values are chosen
    as the 'True' Ih_values, which would be found if unity weights were chosen
    in this example."""
    reflections = flex.reflection_table()
    reflections["intensity"] = flex.double([100.0, 100.0, 80.0, 60.0, 30.0, 40.0, 60.0])
    reflections["inverse_scale_factor"] = flex.double(7, 2.0)
    reflections["variance"] = flex.double([90.0, 100.0, 90.0, 60.0, 30.0, 50.0, 50.0])
    reflections["miller_index"] = flex.miller_index(
        [(1, 0, 0), (0, 0, 1), (-1, 0, 0), (0, 2, 0), (0, 4, 0), (0, 0, -2), (0, 0, 2)]
    )
    reflections.set_flags(flex.bool(7, True), reflections.flags.integrated)
    return reflections


def generate_refl_2():
    """Generate another test reflection table."""
    reflections = flex.reflection_table()
    reflections["intensity"] = flex.double([60.0, 30.0, 10.0, 30.0])
    reflections["variance"] = flex.double([60.0, 30.0, 10.0, 30.0])
    reflections["inverse_scale_factor"] = flex.double(4, 2.0)
    reflections["miller_index"] = flex.miller_index(
        [(1, 0, 0), (0, 4, 0), (10, 0, 0), (0, 4, 0)]
    )
    reflections.set_flags(flex.bool(4, True), reflections.flags.integrated)
    return reflections


def generated_refl_for_splitting_1():
    """ "Make a reflection table."""
    reflections = flex.reflection_table()
    reflections["intensity"] = flex.double([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    reflections["variance"] = flex.double(6, 1.0)
    reflections["inverse_scale_factor"] = flex.double(6, 1.0)
    reflections["miller_index"] = flex.miller_index(
        [(1, 0, 0), (2, 0, 0), (0, 0, 1), (2, 2, 2), (1, 0, 0), (2, 0, 0)]
    )
    reflections["d"] = flex.double([0.8, 2.1, 2.0, 1.4, 1.6, 2.5])
    reflections["partiality"] = flex.double(6, 1.0)
    reflections["xyzobs.px.value"] = flex.vec3_double(
        [
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 5.0),
            (0.0, 0.0, 8.0),
            (0.0, 0.0, 10.0),
            (0.0, 0.0, 12.0),
            (0.0, 0.0, 15.0),
        ]
    )
    reflections["s1"] = flex.vec3_double(
        [
            (0.0, 0.1, 1.0),
            (0.0, 0.1, 1.0),
            (0.0, 0.1, 1.0),
            (0.0, 0.1, 1.0),
            (0.0, 0.1, 1.0),
            (0.0, 0.1, 1.0),
        ]
    )
    reflections.set_flags(flex.bool(6, True), reflections.flags.integrated)
    reflections.set_flags(flex.bool(6, False), reflections.flags.bad_for_scaling)
    return reflections


def generated_refl_for_splitting_2():
    """Make a reflection table."""
    reflections = flex.reflection_table()
    reflections["intensity"] = flex.double([7.0, 8.0, 9.0, 10.0, 11.0])
    reflections["variance"] = flex.double(5, 1.0)
    reflections["inverse_scale_factor"] = flex.double(5, 1.0)
    reflections["miller_index"] = flex.miller_index(
        [(2, 2, 2), (2, 0, 0), (0, 0, 1), (2, 2, 2), (1, 0, 0)]
    )
    reflections["d"] = flex.double([0.8, 2.1, 2.0, 1.4, 1.6])
    reflections["partiality"] = flex.double(5, 1.0)
    reflections["xyzobs.px.value"] = flex.vec3_double(
        [
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 5.0),
            (0.0, 0.0, 8.0),
            (0.0, 0.0, 10.0),
            (0.0, 0.0, 12.0),
        ]
    )
    reflections["s1"] = flex.vec3_double(
        [
            (0.0, 0.1, 1.0),
            (0.0, 0.1, 1.0),
            (0.0, 0.1, 1.0),
            (0.0, 0.1, 1.0),
            (0.0, 0.1, 1.0),
        ]
    )
    reflections.set_flags(flex.bool(5, True), reflections.flags.integrated)
    reflections.set_flags(flex.bool(5, False), reflections.flags.bad_for_scaling)
    return reflections


def test_IhTableblock_onedataset(large_reflection_table, test_sg):
    """Test direct initialisation of Ih_table block"""
    asu_indices = map_indices_to_asu(large_reflection_table["miller_index"], test_sg)
    large_reflection_table["asu_miller_index"] = asu_indices
    n_groups = len(set(asu_indices))
    n_refl = large_reflection_table.size()
    block = IhTableBlock(n_groups=n_groups, n_refl=n_refl)
    group_ids = np.array([0, 1, 0, 2, 3, 4, 4], dtype=np.uint64)

    df = refl_table_to_df(large_reflection_table)

    block.add_data(0, group_ids, df, large_reflection_table["miller_index"])

    # test properties
    assert block.inverse_scale_factors.size
    assert block.variances.size
    assert block.intensities.size
    assert block.weights.size
    assert block.asu_miller_index.size

    # assert list(block.block_selections[0]) == [0, 1, 2, 3, 4, 5, 6]
    assert len(block.block_selections) == 1

    assert block.h_index_matrix.n_rows == n_refl
    assert block.h_index_matrix.n_cols == n_groups
    assert block.h_index_matrix.non_zeroes == n_refl
    assert block.h_index_matrix[0, 0] == 1
    assert block.h_index_matrix[1, 1] == 1
    assert block.h_index_matrix[2, 0] == 1
    assert block.h_index_matrix[3, 2] == 1
    assert block.h_index_matrix[4, 3] == 1
    assert block.h_index_matrix[5, 4] == 1
    assert block.h_index_matrix[6, 4] == 1

    assert block.h_expand_matrix.n_cols == n_refl
    assert block.h_expand_matrix.n_rows == n_groups
    assert block.h_expand_matrix.non_zeroes == n_refl
    assert block.h_expand_matrix[0, 0] == 1
    assert block.h_expand_matrix[1, 1] == 1
    assert block.h_expand_matrix[0, 2] == 1
    assert block.h_expand_matrix[2, 3] == 1
    assert block.h_expand_matrix[3, 4] == 1
    assert block.h_expand_matrix[4, 5] == 1
    assert block.h_expand_matrix[4, 6] == 1

    assert list(block.intensities) == list(large_reflection_table["intensity"])
    assert list(block.weights) == list(1.0 / large_reflection_table["variance"])
    assert "Ih_values" not in block.Ih_table
    block.calc_Ih()
    assert list(block.Ih_values) == pytest.approx(
        [x / 2.0 for x in [90.0, 100.0, 90.0, 60.0, 30.0, 50.0, 50.0]]
    )

    with pytest.raises(AssertionError):
        block.add_data(
            0, group_ids, large_reflection_table, large_reflection_table["miller_index"]
        )

    assert list(block.calc_nh()) == [2, 1, 2, 1, 1, 2, 2]

    # Test update error model
    block.update_weights(mock_error_model())
    assert list(block.weights) == pytest.approx(
        [1.0, 1.0 / 2.0, 1.0 / 3.0, 1.0 / 4.0, 1.0 / 5.0, 1.0 / 6.0, 1.0 / 7.0]
    )

    # Test select method
    new_block = block.select(np.array([True, False, True, False, True, False, False]))

    assert list(new_block.block_selections[0]) == [0, 2, 4]

    assert new_block.h_index_matrix.n_rows == 3
    assert new_block.h_index_matrix.n_cols == 2
    assert new_block.h_index_matrix.non_zeroes == 3
    assert new_block.h_index_matrix[0, 0] == 1
    assert new_block.h_index_matrix[1, 0] == 1
    assert new_block.h_index_matrix[2, 1] == 1

    assert new_block.h_expand_matrix.n_cols == 3
    assert new_block.h_expand_matrix.n_rows == 2
    assert new_block.h_expand_matrix.non_zeroes == 3
    assert new_block.h_expand_matrix[0, 0] == 1
    assert new_block.h_expand_matrix[0, 1] == 1
    assert new_block.h_expand_matrix[1, 2] == 1

    assert list(new_block.Ih_values) == pytest.approx(
        [x / 2.0 for x in [90.0, 90.0, 30.0]]
    )
    with pd.option_context("mode.chained_assignment", None):
        # test setter methods
        new_block.inverse_scale_factors = np.array([1.0, 2.0, 3.0])
        assert list(new_block.inverse_scale_factors) == [1.0, 2.0, 3.0]
        with pytest.raises(AssertionError):
            new_block.inverse_scale_factors = np.array([1.0, 2.0, 3.0, 4.0])
        new_block.weights = np.array([0.1, 0.2, 0.3])
        assert list(new_block.weights) == [0.1, 0.2, 0.3]
        with pytest.raises(AssertionError):
            new_block.weights = np.array([1.0, 2.0, 3.0, 4.0])
        # test selecting initial table on groups, expecting same result as before
        new_block = block.select_on_groups(np.array([True, False, False, True, False]))
        assert list(new_block.block_selections[0]) == [0, 2, 4]

        assert new_block.h_index_matrix.n_rows == 3
        assert new_block.h_index_matrix.n_cols == 2
        assert new_block.h_index_matrix.non_zeroes == 3
        assert new_block.h_index_matrix[0, 0] == 1
        assert new_block.h_index_matrix[1, 0] == 1
        assert new_block.h_index_matrix[2, 1] == 1

        assert new_block.h_expand_matrix.n_cols == 3
        assert new_block.h_expand_matrix.n_rows == 2
        assert new_block.h_expand_matrix.non_zeroes == 3
        assert new_block.h_expand_matrix[0, 0] == 1
        assert new_block.h_expand_matrix[0, 1] == 1
        assert new_block.h_expand_matrix[1, 2] == 1

        assert list(new_block.Ih_values) == pytest.approx(
            [x / 2.0 for x in [90.0, 90.0, 30.0]]
        )


def test_IhTableblock_twodatasets(large_reflection_table, test_sg):
    """Test direct initialisation of Ih_table block. Use the same reflection
    table as previously to make comparisons easier"""
    asu_indices = map_indices_to_asu(large_reflection_table["miller_index"], test_sg)
    large_reflection_table["asu_miller_index"] = asu_indices
    n_groups = len(set(asu_indices))
    n_refl = large_reflection_table.size()

    sel_1 = flex.bool([True, True, False, False, True, True, True])
    dataset_1 = large_reflection_table.select(sel_1)
    dataset_2 = large_reflection_table.select(~sel_1)

    block = IhTableBlock(n_groups=n_groups, n_refl=n_refl, n_datasets=2)

    block.add_data(
        0,
        np.array([0, 1, 3, 4, 4], dtype=np.uint64),
        refl_table_to_df(dataset_1),
        dataset_1["miller_index"],
    )
    block.add_data(
        1,
        np.array([0, 2], dtype=np.uint64),
        refl_table_to_df(dataset_2),
        dataset_2["miller_index"],
    )

    assert list(block.block_selections[0]) == [0, 1, 2, 3, 4]
    assert list(block.block_selections[1]) == [0, 1]
    assert len(block.block_selections) == 2

    assert block.h_index_matrix.n_rows == n_refl
    assert block.h_index_matrix.n_cols == n_groups
    assert block.h_index_matrix.non_zeroes == n_refl
    # first dataset
    assert block.h_index_matrix[0, 0] == 1
    assert block.h_index_matrix[1, 1] == 1
    assert block.h_index_matrix[2, 3] == 1
    assert block.h_index_matrix[3, 4] == 1
    assert block.h_index_matrix[4, 4] == 1
    # then second dataset
    assert block.h_index_matrix[5, 0] == 1
    assert block.h_index_matrix[6, 2] == 1

    assert block.h_expand_matrix.n_cols == n_refl
    assert block.h_expand_matrix.n_rows == n_groups
    assert block.h_expand_matrix.non_zeroes == n_refl
    assert block.h_expand_matrix[0, 0] == 1
    assert block.h_expand_matrix[1, 1] == 1
    assert block.h_expand_matrix[3, 2] == 1
    assert block.h_expand_matrix[4, 3] == 1
    assert block.h_expand_matrix[4, 4] == 1
    assert block.h_expand_matrix[0, 5] == 1
    assert block.h_expand_matrix[2, 6] == 1

    # Test select method
    new_block = block.select(np.array([True, False, True, False, False, True, False]))

    assert list(new_block.block_selections[0]) == [0, 2]
    assert list(new_block.block_selections[1]) == [0]

    assert new_block.h_index_matrix.n_rows == 3
    assert new_block.h_index_matrix.n_cols == 2
    assert new_block.h_index_matrix.non_zeroes == 3
    assert new_block.h_index_matrix[0, 0] == 1
    assert new_block.h_index_matrix[1, 1] == 1
    assert new_block.h_index_matrix[2, 0] == 1

    assert new_block.h_expand_matrix.n_cols == 3
    assert new_block.h_expand_matrix.n_rows == 2
    assert new_block.h_expand_matrix.non_zeroes == 3
    assert new_block.h_expand_matrix[0, 0] == 1
    assert new_block.h_expand_matrix[1, 1] == 1
    assert new_block.h_expand_matrix[0, 2] == 1


def test_IhTable_split_into_blocks(
    large_reflection_table, small_reflection_table, test_sg
):
    """Test that the Ih_table datastructure correctly organises the data
    from two reflection tables into two IhTableBlocks."""

    sel1 = flex.bool(7, True)
    sel1[6] = False
    sel2 = flex.bool(4, True)
    sel2[1] = False

    Ih_table = IhTable(
        reflection_tables=[
            large_reflection_table.select(sel1),
            small_reflection_table.select(sel2),
        ],
        indices_lists=[sel1.iselection(), sel2.iselection()],
        space_group=test_sg,
        nblocks=2,
    )

    assert Ih_table.n_datasets == 2
    assert Ih_table.n_work_blocks == 2
    block_list = Ih_table.Ih_table_blocks
    assert list(block_list[0].asu_miller_index) == [
        (0, 0, 1),
        (0, 0, 2),
        (0, 2, 0),
    ]
    assert list(block_list[1].asu_miller_index) == [
        (0, 4, 0),
        (1, 0, 0),
        (1, 0, 0),
        (0, 4, 0),
        (1, 0, 0),
        (10, 0, 0),
    ]
    assert list(block_list[0].block_selections[0]) == [1, 5, 3]
    assert list(block_list[0].block_selections[1]) == []
    assert list(block_list[1].block_selections[0]) == [4, 0, 2]
    assert list(block_list[1].block_selections[1]) == [3, 0, 2]

    # test the 'get_block_selections_for_dataset' method
    block_sels_0 = Ih_table.get_block_selections_for_dataset(dataset=0)
    assert len(block_sels_0) == 2
    assert list(block_sels_0[0]) == [1, 5, 3]
    assert list(block_sels_0[1]) == [4, 0, 2]
    block_sels_1 = Ih_table.get_block_selections_for_dataset(dataset=1)
    assert len(block_sels_1) == 2
    assert list(block_sels_1[0]) == []
    assert list(block_sels_1[1]) == [3, 0, 2]

    # test the size method
    assert Ih_table.size == 9

    expected_h_idx_matrix = sparse.matrix(4, 3)
    expected_h_idx_matrix[0, 0] = 1
    expected_h_idx_matrix[1, 1] = 1
    expected_h_idx_matrix[2, 2] = 1
    assert block_list[0].h_index_matrix == expected_h_idx_matrix
    expected_h_idx_matrix = sparse.matrix(7, 3)
    expected_h_idx_matrix[0, 0] = 1
    expected_h_idx_matrix[1, 1] = 1
    expected_h_idx_matrix[2, 1] = 1
    expected_h_idx_matrix[3, 0] = 1
    expected_h_idx_matrix[4, 1] = 1
    expected_h_idx_matrix[5, 2] = 1
    assert block_list[1].h_index_matrix == expected_h_idx_matrix

    s = flex.double([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    assert list(block_list[1].sum_in_groups(s, "per_group")) == [5.0, 10.0, 6.0]
    assert list(block_list[1].sum_in_groups(s, "per_refl")) == [
        5.0,
        10.0,
        10.0,
        5.0,
        10.0,
        6.0,
    ]
    with pytest.raises(ValueError):
        block_list[1].sum_in_groups(s, "invalid")

    # check that only dataset 1 intensities are updated
    new_intensities = np.array([60.0, 50.0, 40.0, 30.0, 20.0, 10.0])
    Ih_table.update_data_in_blocks(
        data=new_intensities, dataset_id=0, column="intensity"
    )
    assert list(Ih_table.blocked_data_list[0].intensities) == [50.0, 10.0, 30.0]
    assert list(Ih_table.blocked_data_list[1].intensities) == [
        20.0,
        60.0,
        40.0,
        30.0,
        60.0,
        10.0,
    ]
    # try updating variances
    new_vars = np.array([100.0, 200.0, 300.0, 400.0])
    Ih_table.update_data_in_blocks(data=new_vars, dataset_id=1, column="variance")
    assert list(Ih_table.blocked_data_list[0].variances) == [100.0, 50.0, 60.0]
    assert list(Ih_table.blocked_data_list[1].variances) == [
        30.0,
        90.0,
        90.0,
        400.0,
        100.0,
        300.0,
    ]


def test_IhTable_freework(large_reflection_table, small_reflection_table, test_sg):
    sel1 = flex.bool(7, True)
    sel1[6] = False
    sel2 = flex.bool(4, True)
    sel2[1] = False

    Ih_table = IhTable(
        reflection_tables=[
            large_reflection_table.select(sel1),
            small_reflection_table.select(sel2),
        ],
        indices_lists=[sel1.iselection(), sel2.iselection()],
        space_group=test_sg,
        nblocks=2,
        free_set_percentage=50.0,
    )

    assert len(Ih_table.blocked_data_list) == 3
    assert Ih_table.n_datasets == 2
    assert Ih_table.n_work_blocks == 2
    block_list = Ih_table.Ih_table_blocks

    # two standard blocks
    assert block_list[0].h_index_matrix[0, 0] == 1
    assert block_list[0].h_index_matrix.non_zeroes == 1
    assert block_list[1].h_index_matrix[0, 0] == 1
    assert block_list[1].h_index_matrix[1, 0] == 1
    assert block_list[1].h_index_matrix[2, 1] == 1
    assert block_list[1].h_index_matrix.non_zeroes == 3
    # free set block
    assert block_list[2].h_index_matrix[0, 0] == 1
    assert block_list[2].h_index_matrix[1, 1] == 1
    assert block_list[2].h_index_matrix[2, 2] == 1
    assert block_list[2].h_index_matrix[3, 2] == 1
    assert block_list[2].h_index_matrix[4, 2] == 1
    assert block_list[2].h_index_matrix.non_zeroes == 5

    assert list(block_list[0].block_selections[0]) == [5]
    assert list(block_list[0].block_selections[1]) == []
    assert list(block_list[1].block_selections[0]) == [4]
    assert list(block_list[1].block_selections[1]) == [3, 2]
    assert list(block_list[2].block_selections[0]) == [1, 3, 0, 2]
    assert list(block_list[2].block_selections[1]) == [0]

    # test get_block_selections_for_dataset
    block_sels_0 = Ih_table.get_block_selections_for_dataset(0)
    assert len(block_sels_0) == 3
    assert list(block_sels_0[0]) == [5]
    assert list(block_sels_0[1]) == [4]
    assert list(block_sels_0[2]) == [1, 3, 0, 2]
    block_sels_1 = Ih_table.get_block_selections_for_dataset(1)
    assert len(block_sels_1) == 3
    assert list(block_sels_1[0]) == []
    assert list(block_sels_1[1]) == [3, 2]
    assert list(block_sels_1[2]) == [0]
    with pytest.raises(AssertionError):
        _ = Ih_table.get_block_selections_for_dataset(2)

    Ih_table.calc_Ih()

    # test setting data
    # set scale factors
    new_s_block_2 = np.arange(start=1, stop=6)
    Ih_table.set_inverse_scale_factors(new_s_block_2, 2)
    assert list(Ih_table.Ih_table_blocks[2].inverse_scale_factors) == list(
        new_s_block_2
    )
    # set derivatives
    derivs = Mock()
    Ih_table.set_derivatives(derivs, 0)
    assert Ih_table.Ih_table_blocks[0].derivatives is derivs

    def update_vars_side_effect(*args):
        return np.array([0.5] * len(args[0]))

    # test setting an error model
    em = Mock()
    em.update_variances.side_effect = update_vars_side_effect

    Ih_table.update_weights(em)
    for block in Ih_table.Ih_table_blocks:
        assert list(block.weights) == pytest.approx([2.0] * block.size)

    Ih_table.calc_Ih(1)

    # now test free set with offset
    Ih_table = IhTable(
        reflection_tables=[
            large_reflection_table.select(sel1),
            small_reflection_table.select(sel2),
        ],
        indices_lists=[sel1.iselection(), sel2.iselection()],
        space_group=test_sg,
        nblocks=2,
        free_set_percentage=50.0,
        free_set_offset=1,
    )
    assert len(Ih_table.blocked_data_list) == 3
    assert Ih_table.n_datasets == 2
    assert Ih_table.n_work_blocks == 2
    block_list = Ih_table.Ih_table_blocks

    # two standard blocks
    assert block_list[0].h_index_matrix[0, 0] == 1
    assert block_list[0].h_index_matrix[1, 1] == 1
    assert block_list[0].h_index_matrix.non_zeroes == 2
    assert block_list[1].h_index_matrix[0, 0] == 1
    assert block_list[1].h_index_matrix[1, 0] == 1
    assert block_list[1].h_index_matrix[2, 0] == 1
    assert block_list[1].h_index_matrix.non_zeroes == 3
    # free set block
    assert block_list[2].h_index_matrix[0, 0] == 1
    assert block_list[2].h_index_matrix[1, 1] == 1
    assert block_list[2].h_index_matrix[2, 1] == 1
    assert block_list[2].h_index_matrix[3, 2] == 1
    assert block_list[2].h_index_matrix.non_zeroes == 4

    assert list(block_list[0].block_selections[0]) == [1, 3]
    assert list(block_list[0].block_selections[1]) == []
    assert list(block_list[1].block_selections[0]) == [0, 2]
    assert list(block_list[1].block_selections[1]) == [0]
    assert list(block_list[2].block_selections[0]) == [5, 4]
    assert list(block_list[2].block_selections[1]) == [3, 2]

    Ih_table.calc_Ih()

    # Test the 'as_miller_array' method.
    unit_cell = uctbx.unit_cell((1.0, 1.0, 1.0, 90.0, 90.0, 90.0))

    arr = Ih_table.as_miller_array(unit_cell)
    assert arr.size() == 5
    assert list(arr.indices()) == [
        (0, 0, 1),
        (0, 2, 0),
        (1, 0, 0),
        (1, 0, 0),
        (1, 0, 0),
    ]

    arr = Ih_table.as_miller_array(unit_cell, return_free_set_data=True)
    assert arr.size() == 4
    assert list(arr.indices()) == [(0, 0, 2), (0, 4, 0), (0, 4, 0), (10, 0, 0)]


def test_set_Ih_values_to_target(test_sg):
    """Test the setting of Ih values for targeted scaling."""
    """Generate input for testing joint_Ih_table."""
    Ih_table = IhTable(
        [generated_refl_for_splitting_1(), generated_refl_for_splitting_2()],
        test_sg,
        nblocks=2,
    )
    Ih_table.calc_Ih()
    # First check that values are set up correctly.
    block_list = Ih_table.blocked_data_list
    assert list(block_list[0].Ih_values) == pytest.approx(
        [6.0, 17.0 / 3.0, 17.0 / 3.0, 6.0, 17.0 / 3.0]
    )
    assert list(block_list[1].Ih_values) == pytest.approx(
        [16.0 / 3.0, 16.0 / 3.0, 7.0, 16.0 / 3.0, 7.0, 7.0]
    )

    # set some values in the target
    # change the (2, 0, 0) reflections to (4, 0, 0) to test if they are removed
    # from the blocks
    t1 = generated_refl_for_splitting_1()
    t2 = generated_refl_for_splitting_2()
    t1["miller_index"][1] = (4, 0, 0)
    t1["miller_index"][5] = (4, 0, 0)
    t2["miller_index"][1] = (4, 0, 0)
    target = IhTable([t1, t2], test_sg, nblocks=1)
    target.blocked_data_list[0].Ih_table["Ih_values"] = np.array(
        [0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.1, 0.2, 0.3, 0.4, 0.4]
    )

    for block in block_list:
        block.match_Ih_values_to_target(target)
    assert list(block_list[0].Ih_values) == [0.1, 0.2, 0.2, 0.1, 0.2]
    assert list(block_list[1].Ih_values) == [0.4, 0.4, 0.4]


'''@pytest.mark.xfail(reason='not yet updated code')
def test_apply_iterative_weighting(reflection_table_for_block, test_sg):
  """Test the setting of iterative weights."""

  Ihtableblock = IhTableBlock(reflection_table_for_block[0], test_sg,
    weighting_scheme='GM')

  # After construction, weights should initially be set to inverse variances,
  # to allow the first calculation of Ih (by least-squares approach).
  assert list(Ihtableblock.weights) == [1.0] * 6

  # Now update weights
  Ihtableblock.update_weights()
  assert list(Ihtableblock.weights) != [1.0] * 6 # Check changed, or that
  # new weights are not by chance identical to old weights.
  gIh = Ihtableblock.inverse_scale_factors * Ihtableblock.Ih_values
  t = (Ihtableblock.intensities - gIh) / gIh
  assert list(Ihtableblock.weights) == list(1.0/(1.0 + t**2)**2)


  Ih_table = IhTable([reflection_table_for_block[0]], test_sg,
    weighting_scheme='GM')
  block = Ih_table.blocked_data_list[0]
  assert list(block.weights) == [1.0] * 6
  Ih_table.update_weights()
  gIh = block.inverse_scale_factors * block.Ih_values
  t = (block.intensities - gIh) / gIh
  assert list(block.weights) == list(1.0/(1.0 + t**2)**2)'''
