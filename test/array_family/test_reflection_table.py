from __future__ import absolute_import, division, print_function

import copy
import os
import random

import pytest
from cctbx import sgtbx
from dials.array_family import flex
from dials.util import Sorry
from dxtbx.model import ExperimentList, Experiment, Crystal
from dxtbx.serialize import load


def test_init():
    # test default
    table = flex.reflection_table()
    assert table.is_consistent()
    assert table.nrows() == 0
    assert table.ncols() == 0
    assert table.empty()

    # test with nrows
    table = flex.reflection_table(10)
    assert table.is_consistent()
    assert table.nrows() == 10
    assert table.ncols() == 0
    assert table.empty()

    # test with valid columns
    table = flex.reflection_table(
        [
            ("col1", flex.int(10)),
            ("col2", flex.double(10)),
            ("col3", flex.std_string(10)),
        ]
    )
    assert table.is_consistent()
    assert table.nrows() == 10
    assert table.ncols() == 3
    assert not table.empty()

    # test with invalid columns
    with pytest.raises(RuntimeError):
        _ = flex.reflection_table(
            [
                ("col1", flex.int(10)),
                ("col2", flex.double(20)),
                ("col3", flex.std_string(10)),
            ]
        )


def test_resizing():
    # Create a table with 2 empty columns
    table = flex.reflection_table()
    assert table.empty()
    table["col1"] = flex.int()
    table["col2"] = flex.double()
    assert table.nrows() == 0
    assert table.ncols() == 2
    assert not table.empty()
    assert "col1" in table
    assert "col2" in table
    assert "col3" not in table

    # Create a table with 2 columns and 10 rows
    table = flex.reflection_table()
    table["col1"] = flex.int(10)
    table["col2"] = flex.double(10)
    assert table.nrows() == 10
    assert table.ncols() == 2

    # Add an extra column with the wrong size (throw)
    with pytest.raises(RuntimeError):
        table["col3"] = flex.std_string(20)
    assert table.nrows() == 10
    assert table.ncols() == 2
    assert table.is_consistent()
    assert len(table["col1"]) == 10
    assert len(table["col2"]) == 10
    assert len(table) == table.size()

    # Resize the table (should resize all columns)
    table.resize(50)
    assert table.nrows() == 50
    assert table.ncols() == 2
    assert table.is_consistent()
    assert len(table["col1"]) == 50
    assert len(table["col2"]) == 50

    # Make the table inconsistent
    table["col1"].resize(40)
    assert not table.is_consistent()
    with pytest.raises(Exception):
        table.nrows()
    assert table.ncols() == 2

    # Clear the table
    table.clear()
    assert table.is_consistent()
    assert table.empty()
    assert table.nrows() == 0
    assert table.ncols() == 0


def test_delete():
    # Test del item
    table = flex.reflection_table()
    table["col1"] = flex.int([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    table["col2"] = flex.int([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    table["col3"] = flex.int([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    del table["col3"]
    assert table.is_consistent()
    assert table.nrows() == 10
    assert table.ncols() == 2
    assert "col3" not in table

    # Test del row
    del table[5]
    assert table.is_consistent()
    assert table.nrows() == 9
    assert table.ncols() == 2
    assert all(a == b for a, b in zip(list(table["col1"]), [0, 1, 2, 3, 4, 6, 7, 8, 9]))

    # Test del slice
    del table[0:10:2]
    assert table.is_consistent()
    assert table.nrows() == 4
    assert table.ncols() == 2
    assert all(a == b for a, b in zip(list(table["col1"]), [1, 3, 6, 8]))

    # Test del slice
    del table[:]
    assert table.is_consistent()
    assert table.nrows() == 0
    assert table.ncols() == 2


def test_row_operations():
    # The columns as lists
    c1 = list(range(10))
    c2 = list(range(10))
    c3 = ["a", "b", "c", "d", "e", "f", "g", "i", "j", "k"]

    # Create a table with some elements
    table = flex.reflection_table()
    table["col1"] = flex.int(c1)
    table["col2"] = flex.double(c2)
    table["col3"] = flex.std_string(c3)

    # Extend the table
    table.extend(table)
    c1 = c1 * 2
    c2 = c2 * 2
    c3 = c3 * 2
    assert table.nrows() == 20
    assert table.ncols() == 3
    assert table.is_consistent()
    assert all(a == b for a, b in zip(table["col1"], c1))
    assert all(a == b for a, b in zip(table["col2"], c2))
    assert all(a == b for a, b in zip(table["col3"], c3))

    # Append some rows to the table
    row = {"col1": 10}
    c1 = c1 + [10]
    c2 = c2 + [0]
    c3 = c3 + [""]
    table.append(row)
    assert table.nrows() == 21
    assert table.ncols() == 3
    assert table.is_consistent()
    assert all(a == b for a, b in zip(table["col1"], c1))
    assert all(a == b for a, b in zip(table["col2"], c2))
    assert all(a == b for a, b in zip(table["col3"], c3))

    row = {"col2": 11}
    c1 = c1 + [0]
    c2 = c2 + [11]
    c3 = c3 + [""]
    table.append(row)
    assert table.nrows() == 22
    assert table.ncols() == 3
    assert table.is_consistent()
    assert all(a == b for a, b in zip(table["col1"], c1))
    assert all(a == b for a, b in zip(table["col2"], c2))
    assert all(a == b for a, b in zip(table["col3"], c3))

    row = {"col1": 12, "col2": 12, "col3": "l"}
    c1 = c1 + [12]
    c2 = c2 + [12]
    c3 = c3 + ["l"]
    table.append(row)
    assert table.nrows() == 23
    assert table.ncols() == 3
    assert table.is_consistent()
    assert all(a == b for a, b in zip(table["col1"], c1))
    assert all(a == b for a, b in zip(table["col2"], c2))
    assert all(a == b for a, b in zip(table["col3"], c3))

    # Try inserting some rows
    row = {"col1": -1}
    c1.insert(5, -1)
    c2.insert(5, 0)
    c3.insert(5, "")
    table.insert(5, row)
    assert table.nrows() == 24
    assert table.ncols() == 3
    assert table.is_consistent()
    assert all(a == b for a, b in zip(table["col1"], c1))
    assert all(a == b for a, b in zip(table["col2"], c2))
    assert all(a == b for a, b in zip(table["col3"], c3))

    row = {"col1": -2, "col2": -3, "col3": "abc"}
    c1.insert(2, -2)
    c2.insert(2, -3)
    c3.insert(2, "abc")
    table.insert(2, row)
    assert table.nrows() == 25
    assert table.ncols() == 3
    assert table.is_consistent()
    assert all(a == b for a, b in zip(table["col1"], c1))
    assert all(a == b for a, b in zip(table["col2"], c2))
    assert all(a == b for a, b in zip(table["col3"], c3))

    # Try iterating through table rows
    for i in range(table.nrows()):
        row = table[i]
        assert row["col1"] == c1[i]
        assert row["col2"] == c2[i]
        assert row["col3"] == c3[i]

    # Trying setting some rows
    row = {"col1": 100}
    table[2] = row
    assert table[2]["col1"] == 100
    assert table[2]["col2"] == c2[2]
    assert table[2]["col3"] == c3[2]

    row = {"col1": 1000, "col2": 2000, "col3": "hello"}
    table[10] = row
    assert table[10]["col1"] == 1000
    assert table[10]["col2"] == 2000
    assert table[10]["col3"] == "hello"


def test_iteration():
    # The columns as lists
    c1 = list(range(10))
    c2 = list(range(10))
    c3 = ["a", "b", "c", "d", "e", "f", "g", "i", "j", "k"]

    # Create a table with some elements
    table = flex.reflection_table()
    table["col1"] = flex.int(c1)
    table["col2"] = flex.double(c2)
    table["col3"] = flex.std_string(c3)

    # Try iterating keys
    k = []
    for key in table.keys():
        k.append(key)
    assert len(k) == 3
    assert k.count("col1") == 1
    assert k.count("col2") == 1
    assert k.count("col3") == 1

    # Try iterating columns
    k = []
    c = []
    for key, col in table.cols():
        k.append(key)
        c.append(col)
    assert len(k) == 3
    assert k.count("col1") == 1
    assert k.count("col2") == 1
    assert k.count("col3") == 1

    # Try iterating rows
    for row1, row2 in zip(table.rows(), zip(c1, c2, c3)):
        assert row1["col1"] == row2[0]
        assert row1["col2"] == row2[1]
        assert row1["col3"] == row2[2]


def test_slicing():
    # The columns as lists
    c1 = list(range(10))
    c2 = list(range(10))
    c3 = ["a", "b", "c", "d", "e", "f", "g", "i", "j", "k"]

    # Create a table with some elements
    table = flex.reflection_table()
    table["col1"] = flex.int(c1)
    table["col2"] = flex.double(c2)
    table["col3"] = flex.std_string(c3)

    # Try forward slicing
    new_table = table[2:7:2]
    assert new_table.ncols() == 3
    assert new_table.nrows() == 3
    assert new_table.is_consistent()
    c11 = c1[2:7:2]
    c22 = c2[2:7:2]
    c33 = c3[2:7:2]
    assert all(a == b for a, b in zip(new_table["col1"], c11))
    assert all(a == b for a, b in zip(new_table["col2"], c22))
    assert all(a == b for a, b in zip(new_table["col3"], c33))

    # Try backward slicing
    new_table = table[7:2:-2]
    assert new_table.ncols() == 3
    assert new_table.nrows() == 3
    assert new_table.is_consistent()
    c11 = c1[7:2:-2]
    c22 = c2[7:2:-2]
    c33 = c3[7:2:-2]
    assert all(a == b for a, b in zip(new_table["col1"], c11))
    assert all(a == b for a, b in zip(new_table["col2"], c22))
    assert all(a == b for a, b in zip(new_table["col3"], c33))

    # Try setting forward slicing
    table[2:7:2] = new_table
    assert table.ncols() == 3
    assert table.nrows() == 10
    assert table.is_consistent()
    c1[2:7:2] = c11
    c2[2:7:2] = c22
    c3[2:7:2] = c33
    assert all(a == b for a, b in zip(table["col1"], c1))
    assert all(a == b for a, b in zip(table["col2"], c2))
    assert all(a == b for a, b in zip(table["col3"], c3))

    # Try setting backward slicing
    table[7:2:-2] = new_table
    assert table.ncols() == 3
    assert table.nrows() == 10
    assert table.is_consistent()
    c1[7:2:-2] = c11
    c2[7:2:-2] = c22
    c3[7:2:-2] = c33
    assert all(a == b for a, b in zip(table["col1"], c1))
    assert all(a == b for a, b in zip(table["col2"], c2))
    assert all(a == b for a, b in zip(table["col3"], c3))


def test_updating():
    # The columns as lists
    c1 = list(range(10))
    c2 = list(range(10))
    c3 = ["a", "b", "c", "d", "e", "f", "g", "i", "j", "k"]

    # Create a table with some elements
    table0 = flex.reflection_table()
    table1 = flex.reflection_table()
    table2 = flex.reflection_table()
    table1["col1"] = flex.int(c1)
    table1["col2"] = flex.double(c2)
    table2["col3"] = flex.std_string(c3)

    # Update from zero columns
    table0.update(table1)
    assert table0.is_consistent()
    assert table0.nrows() == 10
    assert table0.ncols() == 2

    # Update table1 with table2 columns
    table1.update(table2)
    assert table1.is_consistent()
    assert table1.nrows() == 10
    assert table1.ncols() == 3
    assert table2.is_consistent()
    assert table2.nrows() == 10
    assert table2.ncols() == 1

    # Update trable1 with invalid table
    c3 = ["a", "b", "c"]

    # Create a table with some elements
    table2 = flex.reflection_table()
    table2["col3"] = flex.std_string(c3)
    with pytest.raises(RuntimeError):
        table1.update(table2)

    assert table1.is_consistent()
    assert table1.nrows() == 10
    assert table1.ncols() == 3
    assert table2.is_consistent()
    assert table2.nrows() == 3
    assert table2.ncols() == 1


def test_select():
    # The columns as lists
    c1 = list(range(10))
    c2 = list(range(10))
    c3 = ["a", "b", "c", "d", "e", "f", "g", "i", "j", "k"]

    # Create a table with some elements
    table = flex.reflection_table()
    table["col1"] = flex.int(c1)
    table["col2"] = flex.double(c2)
    table["col3"] = flex.std_string(c3)

    # Select some columns
    new_table = table.select(("col1", "col2"))
    assert new_table.nrows() == 10
    assert new_table.ncols() == 2
    assert all(a == b for a, b in zip(new_table["col1"], c1))
    assert all(a == b for a, b in zip(new_table["col2"], c2))

    # Select some columns
    new_table = table.select(flex.std_string(["col1", "col2"]))
    assert new_table.nrows() == 10
    assert new_table.ncols() == 2
    assert all(a == b for a, b in zip(new_table["col1"], c1))
    assert all(a == b for a, b in zip(new_table["col2"], c2))

    # Select some rows
    index = flex.size_t([0, 1, 5, 8, 9])
    cc1 = [c1[i] for i in index]
    cc2 = [c2[i] for i in index]
    cc3 = [c3[i] for i in index]
    new_table = table.select(index)
    assert new_table.nrows() == 5
    assert new_table.ncols() == 3
    assert all(a == b for a, b in zip(new_table["col1"], cc1))
    assert all(a == b for a, b in zip(new_table["col2"], cc2))
    assert all(a == b for a, b in zip(new_table["col3"], cc3))

    # Select some rows
    index = flex.bool([True, True, False, False, False, True, False, False, True, True])
    new_table = table.select(index)
    assert new_table.nrows() == 5
    assert new_table.ncols() == 3
    assert all(a == b for a, b in zip(new_table["col1"], cc1))
    assert all(a == b for a, b in zip(new_table["col2"], cc2))
    assert all(a == b for a, b in zip(new_table["col3"], cc3))


def test_set_selected():
    # The columns as lists
    c1 = list(range(10))
    c2 = list(range(10))
    c3 = ["a", "b", "c", "d", "e", "f", "g", "i", "j", "k"]

    # Create a table with some elements
    table1 = flex.reflection_table()
    table2 = flex.reflection_table()
    table1["col1"] = flex.int(c1)
    table2["col2"] = flex.double(c2)
    table2["col3"] = flex.std_string(c3)

    # Set selected columns
    table1.set_selected(("col3", "col2"), table2)
    assert table1.nrows() == 10
    assert table1.ncols() == 3
    assert all(a == b for a, b in zip(table1["col1"], c1))
    assert all(a == b for a, b in zip(table1["col2"], c2))
    assert all(a == b for a, b in zip(table1["col3"], c3))

    # Set selected columns
    table1 = flex.reflection_table()
    table1["col1"] = flex.int(c1)
    table1.set_selected(flex.std_string(["col3", "col2"]), table2)
    assert table1.nrows() == 10
    assert table1.ncols() == 3
    assert all(a == b for a, b in zip(table1["col1"], c1))
    assert all(a == b for a, b in zip(table1["col2"], c2))
    assert all(a == b for a, b in zip(table1["col3"], c3))

    cc1 = list(range(10, 15))
    cc2 = list(range(10, 15))
    cc3 = ["l", "m", "n", "o", "p"]

    # Set selected rows
    table2 = flex.reflection_table()
    table2["col1"] = flex.int(cc1)
    table2["col2"] = flex.double(cc2)
    table2["col3"] = flex.std_string(cc3)

    index = flex.size_t([0, 1, 5, 8, 9])
    ccc1 = copy.deepcopy(c1)
    ccc2 = copy.deepcopy(c2)
    ccc3 = copy.deepcopy(c3)
    for j, i in enumerate(index):
        ccc1[i] = cc1[j]
        ccc2[i] = cc2[j]
        ccc3[i] = cc3[j]
    table1.set_selected(index, table2)
    assert all(a == b for a, b in zip(table1["col1"], ccc1))
    assert all(a == b for a, b in zip(table1["col2"], ccc2))
    assert all(a == b for a, b in zip(table1["col3"], ccc3))

    # Set selected rows
    table2 = flex.reflection_table()
    table2["col1"] = flex.int(cc1)
    table2["col2"] = flex.double(cc2)
    table2["col3"] = flex.std_string(cc3)

    table1.set_selected(index, table2)
    assert all(a == b for a, b in zip(table1["col1"], ccc1))
    assert all(a == b for a, b in zip(table1["col2"], ccc2))
    assert all(a == b for a, b in zip(table1["col3"], ccc3))


def test_del_selected():
    # The columns as lists
    c1 = list(range(10))
    c2 = list(range(10))
    c3 = ["a", "b", "c", "d", "e", "f", "g", "i", "j", "k"]

    # Create a table with some elements
    table1 = flex.reflection_table()
    table1["col1"] = flex.int(c1)
    table1["col2"] = flex.double(c2)
    table1["col3"] = flex.std_string(c3)

    # Del selected columns
    table1.del_selected(("col3", "col2"))
    assert table1.nrows() == 10
    assert table1.ncols() == 1
    assert "col1" in table1
    assert "col2" not in table1
    assert "col3" not in table1
    assert all(a == b for a, b in zip(table1["col1"], c1))

    # Del selected columns
    table1 = flex.reflection_table()
    table1["col1"] = flex.int(c1)
    table1["col2"] = flex.double(c2)
    table1["col3"] = flex.std_string(c3)
    table1.del_selected(flex.std_string(["col3", "col2"]))
    assert table1.nrows() == 10
    assert table1.ncols() == 1
    assert "col1" in table1
    assert "col2" not in table1
    assert "col3" not in table1
    assert all(a == b for a, b in zip(table1["col1"], c1))

    # Del selected rows
    table1 = flex.reflection_table()
    table1["col1"] = flex.int(c1)
    table1["col2"] = flex.double(c2)
    table1["col3"] = flex.std_string(c3)

    index = flex.size_t([0, 1, 5, 8, 9])
    index2 = list(range(10))
    for i in index:
        index2.remove(i)
    ccc1 = [c1[i] for i in index2]
    ccc2 = [c2[i] for i in index2]
    ccc3 = [c3[i] for i in index2]
    table1.del_selected(index)
    assert table1.nrows() == len(ccc1)
    assert all(a == b for a, b in zip(table1["col1"], ccc1))
    assert all(a == b for a, b in zip(table1["col2"], ccc2))
    assert all(a == b for a, b in zip(table1["col3"], ccc3))

    # Del selected rows
    table1 = flex.reflection_table()
    table1["col1"] = flex.int(c1)
    table1["col2"] = flex.double(c2)
    table1["col3"] = flex.std_string(c3)

    table1.del_selected(index)
    assert table1.nrows() == len(ccc1)
    assert all(a == b for a, b in zip(table1["col1"], ccc1))
    assert all(a == b for a, b in zip(table1["col2"], ccc2))
    assert all(a == b for a, b in zip(table1["col3"], ccc3))


def test_sort():
    table = flex.reflection_table()
    table["a"] = flex.int([2, 4, 3, 1, 5])
    table["b"] = flex.vec2_double([(3, 2), (3, 1), (1, 3), (4, 5), (4, 3)])
    table["c"] = flex.miller_index(
        [(3, 2, 1), (3, 1, 1), (2, 4, 2), (2, 1, 1), (1, 1, 1)]
    )

    table.sort("a")
    assert list(table["a"]) == [1, 2, 3, 4, 5]

    table.sort("b")
    assert list(table["b"]) == [(1, 3), (3, 1), (3, 2), (4, 3), (4, 5)]

    table.sort("c")
    assert list(table["c"]) == [(1, 1, 1), (2, 1, 1), (2, 4, 2), (3, 1, 1), (3, 2, 1)]

    table.sort("c", order=(1, 2, 0))
    assert list(table["c"]) == [(1, 1, 1), (2, 1, 1), (3, 1, 1), (3, 2, 1), (2, 4, 2)]


def test_flags():
    # Create a table with flags all 0
    table = flex.reflection_table()
    table["flags"] = flex.size_t(5, 0)

    # Get all the flags
    f1 = table.get_flags(table.flags.predicted)
    assert f1.count(True) == 0

    # Set some flags
    mask = flex.bool([True, True, False, False, True])
    table.set_flags(mask, table.flags.predicted)
    f1 = table.get_flags(table.flags.predicted)
    assert f1.count(True) == 3
    assert all(f11 == f22 for f11, f22 in zip(f1, mask))
    f2 = table.get_flags(table.flags.predicted | table.flags.observed)
    assert f2.count(True) == 0

    # Unset the flags
    mask = flex.bool(5, True)
    table.unset_flags(mask, table.flags.predicted | table.flags.observed)
    f1 = table.get_flags(table.flags.predicted)
    assert f1.count(True) == 0
    flags = table["flags"]
    assert all(f == 0 for f in flags)

    # Set multiple flags
    mask = flex.bool([True, True, False, False, True])
    table.set_flags(mask, table.flags.predicted | table.flags.observed)
    f1 = table.get_flags(table.flags.predicted)
    f2 = table.get_flags(table.flags.observed)
    assert f1.count(True) == 3
    assert f2.count(True) == 3
    mask = flex.bool([False, True, True, True, False])
    table.set_flags(mask, table.flags.integrated)
    f1 = table.get_flags(table.flags.predicted)
    f2 = table.get_flags(table.flags.observed)
    f3 = table.get_flags(table.flags.integrated)
    f4 = table.get_flags(table.flags.integrated | table.flags.predicted)
    assert f1.count(True) == 3
    assert f2.count(True) == 3
    assert f3.count(True) == 3
    assert f4.count(True) == 1

    # Get where any are set
    f1 = table.get_flags(table.flags.predicted, all=False)
    f2 = table.get_flags(table.flags.observed, all=False)
    f3 = table.get_flags(table.flags.integrated, all=False)
    f4 = table.get_flags(table.flags.integrated | table.flags.predicted, all=False)
    assert f1.count(True) == 3
    assert f2.count(True) == 3
    assert f3.count(True) == 3
    assert f4.count(True) == 5


def test_serialize():
    # The columns as lists
    c1 = list(range(10))
    c2 = list(range(10))
    c3 = ["a", "b", "c", "d", "e", "f", "g", "i", "j", "k"]

    # Create a table with some elements
    table = flex.reflection_table()
    table["col1"] = flex.int(c1)
    table["col2"] = flex.double(c2)
    table["col3"] = flex.std_string(c3)

    # Pickle, then unpickle
    import six.moves.cPickle as pickle

    obj = pickle.dumps(table)
    new_table = pickle.loads(obj)
    assert new_table.is_consistent()
    assert new_table.nrows() == 10
    assert new_table.ncols() == 3
    assert all(a == b for a, b in zip(new_table["col1"], c1))
    assert all(a == b for a, b in zip(new_table["col2"], c2))
    assert all(a == b for a, b in zip(new_table["col3"], c3))


def test_copy():
    # Create a table
    table = flex.reflection_table([("col1", flex.int(range(10)))])

    # Make a shallow copy of the table
    shallow = copy.copy(table)
    shallow["col2"] = flex.double(range(10))
    assert table.ncols() == 2
    assert table.is_consistent()

    # Make a deep copy of the table
    deep = copy.deepcopy(table)
    deep["col3"] = flex.std_string(10)
    assert table.ncols() == 2
    assert deep.ncols() == 3
    assert table.is_consistent()
    assert deep.is_consistent()

    table2 = table.copy()
    table2["col3"] = flex.std_string(10)
    assert table.ncols() == 2
    assert table2.ncols() == 3
    assert table.is_consistent()
    assert table2.is_consistent()


def test_extract_shoeboxes():
    from dials.algorithms.shoebox import MaskCode

    random.seed(0)

    reflections = flex.reflection_table()
    reflections["panel"] = flex.size_t()
    reflections["bbox"] = flex.int6()

    npanels = 2
    width = 1000
    height = 1000
    frame0 = 10
    frame1 = 100
    nrefl = 1000

    for i in range(nrefl):
        xs = random.randint(5, 10)
        ys = random.randint(5, 10)
        x0 = random.randint(-xs + 1, width - 1)
        y0 = random.randint(-ys + 1, height - 1)
        z0 = random.randint(frame0, frame1 - 1)
        x1 = x0 + xs
        y1 = y0 + ys
        z1 = min([z0 + random.randint(1, 10), frame1])
        assert x1 > x0
        assert y1 > y0
        assert z1 > z0
        assert z0 >= frame0 and z1 <= frame1
        bbox = (x0, x1, y0, y1, z0, z1)
        reflections.append({"panel": random.randint(0, 1), "bbox": bbox})

    reflections["shoebox"] = flex.shoebox(reflections["panel"], reflections["bbox"])
    reflections["shoebox"].allocate()

    class FakeImageSet(object):
        def __init__(self):
            self.data = flex.int(range(height * width))
            self.data.reshape(flex.grid(height, width))

        def get_array_range(self):
            return (frame0, frame1)

        def get_detector(self):
            class FakeDetector(object):
                def __len__(self):
                    return npanels

                def __getitem__(self, index):
                    class FakePanel(object):
                        def get_trusted_range(self):
                            return (-1, 1000000)

                    return FakePanel()

            return FakeDetector()

        def __len__(self):
            return frame1 - frame0

        def __getitem__(self, index):
            f = frame0 + index
            return (self.data + f * 1, self.data + f * 2)

        def get_corrected_data(self, index):
            f = frame0 + index
            return (self.data + f * 1, self.data + f * 2)

        def get_mask(self, index):
            image = self.get_corrected_data(index)
            return tuple(im >= 0 for im in image)

    imageset = FakeImageSet()

    reflections.extract_shoeboxes(imageset)

    for i in range(len(reflections)):
        sbox = reflections[i]["shoebox"]
        assert sbox.is_consistent()
        mask = sbox.mask
        data = sbox.data
        bbox = sbox.bbox
        panel = sbox.panel
        x0, x1, y0, y1, z0, z1 = bbox
        for z in range(z1 - z0):
            for y in range(y1 - y0):
                for x in range(x1 - x0):
                    v1 = data[z, y, x]
                    m1 = mask[z, y, x]
                    if (
                        x0 + x >= 0
                        and y0 + y >= 0
                        and x0 + x < width
                        and y0 + y < height
                    ):
                        v2 = imageset.data[y + y0, x + x0] + (z + z0) * (panel + 1)
                        m2 = MaskCode.Valid
                        assert v1 == v2
                        assert m1 == m2
                    else:
                        assert v1 == 0
                        assert m1 == 0


def test_split_by_experiment_id():
    r = flex.reflection_table()
    r["id"] = flex.int()
    for i in range(100):
        r.append({"id": 0})
        r.append({"id": 1})
        r.append({"id": 2})
        r.append({"id": 3})
        r.append({"id": 5})
    result = r.split_by_experiment_id()
    assert len(result) == 5
    for res, exp in zip(result, [0, 1, 2, 3, 5]):
        assert len(res) == 100
        assert res["id"].count(exp) == 100

    # test the same but with experiment_identifiers() set - keep separate as
    # function must work with and without experiment_identifiers() set
    r.experiment_identifiers()[0] = "0"
    r.experiment_identifiers()[1] = "1"
    r.experiment_identifiers()[2] = "2"
    r.experiment_identifiers()[3] = "3"
    r.experiment_identifiers()[5] = "5"
    result = r.split_by_experiment_id()
    assert len(result) == 5
    for res, exp in zip(result, [0, 1, 2, 3, 5]):
        assert len(res) == 100
        assert res["id"].count(exp) == 100
        assert list(res.experiment_identifiers().keys()) == [exp]
        assert list(res.experiment_identifiers().values()) == [str(exp)]


def test_split_indices_by_experiment_id():
    r = flex.reflection_table()
    r["id"] = flex.int()
    for i in range(100):
        r.append({"id": 0})
        r.append({"id": 1})
        r.append({"id": 2})
        r.append({"id": 3})
        r.append({"id": 5})
    index_list = r.split_indices_by_experiment_id(6)
    assert len(index_list) == 6
    for index, exp, num in zip(
        index_list, [0, 1, 2, 3, 4, 5], [100, 100, 100, 100, 0, 100]
    ):
        assert len(index) == num
        assert r.select(index)["id"].count(exp) == num


def test_split_partials():
    r = flex.reflection_table()
    r["value1"] = flex.double()
    r["value2"] = flex.int()
    r["value3"] = flex.double()
    r["bbox"] = flex.int6()
    expected = []
    for i in range(100):
        x0 = random.randint(0, 100)
        x1 = x0 + random.randint(1, 10)
        y0 = random.randint(0, 100)
        y1 = y0 + random.randint(1, 10)
        z0 = random.randint(0, 100)
        z1 = z0 + random.randint(1, 10)
        v1 = random.uniform(0, 100)
        v2 = random.randint(0, 100)
        v3 = random.uniform(0, 100)
        r.append(
            {"value1": v1, "value2": v2, "value3": v3, "bbox": (x0, x1, y0, y1, z0, z1)}
        )
        for z in range(z0, z1):
            expected.append(
                {
                    "value1": v1,
                    "value2": v2,
                    "value3": v3,
                    "bbox": (x0, x1, y0, y1, z, z + 1),
                    "partial_id": i,
                }
            )

    r.split_partials()
    assert len(r) == len(expected)
    EPS = 1e-7
    for r1, r2 in zip(r, expected):
        assert abs(r1["value1"] - r2["value1"]) < EPS
        assert r1["value2"] == r2["value2"]
        assert abs(r1["value3"] - r2["value3"]) < EPS
        assert r1["bbox"] == r2["bbox"]
        assert r1["partial_id"] == r2["partial_id"]


def test_split_partials_with_shoebox():
    from dials.model.data import Shoebox

    r = flex.reflection_table()
    r["value1"] = flex.double()
    r["value2"] = flex.int()
    r["value3"] = flex.double()
    r["bbox"] = flex.int6()
    r["panel"] = flex.size_t()
    r["shoebox"] = flex.shoebox()
    expected = []
    for i in range(100):
        x0 = random.randint(0, 100)
        x1 = x0 + random.randint(1, 10)
        y0 = random.randint(0, 100)
        y1 = y0 + random.randint(1, 10)
        z0 = random.randint(0, 100)
        z1 = z0 + random.randint(1, 10)
        v1 = random.uniform(0, 100)
        v2 = random.randint(0, 100)
        v3 = random.uniform(0, 100)
        sbox = Shoebox(0, (x0, x1, y0, y1, z0, z1))
        sbox.allocate()
        assert sbox.is_consistent()
        w = x1 - x0
        h = y1 - y0
        for z in range(z0, z1):
            for y in range(y0, y1):
                for x in range(x0, x1):
                    sbox.data[z - z0, y - y0, x - x0] = x + y * w + z * w * h
        r.append(
            {
                "value1": v1,
                "value2": v2,
                "value3": v3,
                "bbox": (x0, x1, y0, y1, z0, z1),
                "panel": 0,
                "shoebox": sbox,
            }
        )
        for z in range(z0, z1):
            sbox = Shoebox(0, (x0, x1, y0, y1, z, z + 1))
            sbox.allocate()
            assert sbox.is_consistent()
            w = x1 - x0
            h = y1 - y0
            for y in range(y0, y1):
                for x in range(x0, x1):
                    sbox.data[0, y - y0, x - x0] = x + y * w + z * w * h
            expected.append(
                {
                    "value1": v1,
                    "value2": v2,
                    "value3": v3,
                    "bbox": (x0, x1, y0, y1, z, z + 1),
                    "partial_id": i,
                    "panel": 0,
                    "shoebox": sbox,
                }
            )

    r.split_partials_with_shoebox()
    assert len(r) == len(expected)
    EPS = 1e-7
    for r1, r2 in zip(r, expected):
        assert abs(r1["value1"] - r2["value1"]) < EPS
        assert r1["value2"] == r2["value2"]
        assert abs(r1["value3"] - r2["value3"]) < EPS
        assert r1["bbox"] == r2["bbox"]
        assert r1["partial_id"] == r2["partial_id"]
        assert r1["panel"] == r2["panel"]
        assert (
            r1["shoebox"]
            .data.as_double()
            .as_1d()
            .all_approx_equal(r2["shoebox"].data.as_double().as_1d())
        )


def test_find_overlapping():
    N = 10000
    r = flex.reflection_table(N)
    r["bbox"] = flex.int6(N)
    r["panel"] = flex.size_t(N)
    r["id"] = flex.int(N)
    r["imageset_id"] = flex.int(N)
    for i in range(N):
        x0 = random.randint(0, 100)
        x1 = random.randint(1, 10) + x0
        y0 = random.randint(0, 100)
        y1 = random.randint(1, 10) + y0
        z0 = random.randint(0, 100)
        z1 = random.randint(1, 10) + z0
        panel = random.randint(0, 2)
        pid = random.randint(0, 2)
        r["bbox"][i] = (x0, x1, y0, y1, z0, z1)
        r["panel"][i] = panel
        r["id"][i] = pid
        r["imageset_id"][i] = pid

    def is_overlap(b0, b1, border):
        b0 = (
            b0[0] - border,
            b0[1] + border,
            b0[2] - border,
            b0[3] + border,
            b0[4] - border,
            b0[5] + border,
        )
        b1 = (
            b1[0] - border,
            b1[1] + border,
            b1[2] - border,
            b1[3] + border,
            b1[4] - border,
            b1[5] + border,
        )
        if not (
            b1[0] > b0[1]
            or b1[1] < b0[0]
            or b1[2] > b0[3]
            or b1[3] < b0[2]
            or b1[4] > b0[5]
            or b1[5] < b0[4]
        ):
            return True
        return False

    for i in [0, 2, 5]:
        overlaps = r.find_overlaps(border=i)
        for item in overlaps.edges():
            i0 = overlaps.source(item)
            i1 = overlaps.target(item)
            r0 = r[i0]
            r1 = r[i1]
            p0 = r0["panel"]
            p1 = r1["panel"]
            b0 = r0["bbox"]
            b1 = r1["bbox"]
            j0 = r0["imageset_id"]
            j1 = r1["imageset_id"]
            assert j0 == j1
            assert p0 == p1
            assert is_overlap(b0, b1, i)


def test_to_from_msgpack(tmpdir):
    from dials.model.data import Shoebox

    def gen_shoebox():
        shoebox = Shoebox(0, (0, 4, 0, 3, 0, 1))
        shoebox.allocate()
        for k in range(1):
            for j in range(3):
                for i in range(4):
                    shoebox.data[k, j, i] = i + j + k
                    shoebox.mask[k, j, i] = i % 2
                    shoebox.background[k, j, i] = i * j
        return shoebox

    def compare(a, b):
        assert a.panel == b.panel
        assert a.bbox == b.bbox
        for aa, bb in zip(a.data, b.data):
            if abs(aa - bb) > 1e-9:
                return False
        for aa, bb in zip(a.background, b.background):
            if abs(aa - bb) > 1e-9:
                return False
        for aa, bb in zip(a.mask, b.mask):
            if aa != bb:
                return False
        return True

    # The columns as lists
    c1 = list(range(10))
    c2 = list(range(10))
    c3 = ["a", "b", "c", "d", "e", "f", "g", "i", "j", "k"]
    c4 = [True, False, True, False, True] * 2
    c5 = list(range(10))
    c6 = [(i + 1, i + 2) for i in range(10)]
    c7 = [(i + 1, i + 2, i + 3) for i in range(10)]
    c8 = [tuple(i + j for j in range(9)) for i in range(10)]
    c9 = [tuple(i + j for j in range(6)) for i in range(10)]
    c10 = [(i + 1, i + 2, i + 3) for i in range(10)]
    c11 = [gen_shoebox() for i in range(10)]

    # Create a table with some elements
    table = flex.reflection_table()
    table["col1"] = flex.int(c1)
    table["col2"] = flex.double(c2)
    table["col3"] = flex.std_string(c3)
    table["col4"] = flex.bool(c4)
    table["col5"] = flex.size_t(c5)
    table["col6"] = flex.vec2_double(c6)
    table["col7"] = flex.vec3_double(c7)
    table["col8"] = flex.mat3_double(c8)
    table["col9"] = flex.int6(c9)
    table["col10"] = flex.miller_index(c10)
    table["col11"] = flex.shoebox(c11)

    obj = table.as_msgpack()
    new_table = flex.reflection_table.from_msgpack(obj)
    assert new_table.is_consistent()
    assert new_table.nrows() == 10
    assert new_table.ncols() == 11
    assert all(tuple(a == b for a, b in zip(new_table["col1"], c1)))
    assert all(tuple(a == b for a, b in zip(new_table["col2"], c2)))
    assert all(tuple(a == b for a, b in zip(new_table["col3"], c3)))
    assert all(tuple(a == b for a, b in zip(new_table["col4"], c4)))
    assert all(tuple(a == b for a, b in zip(new_table["col5"], c5)))
    assert all(tuple(a == b for a, b in zip(new_table["col6"], c6)))
    assert all(tuple(a == b for a, b in zip(new_table["col7"], c7)))
    assert all(tuple(a == b for a, b in zip(new_table["col8"], c8)))
    assert all(tuple(a == b for a, b in zip(new_table["col9"], c9)))
    assert all(tuple(a == b for a, b in zip(new_table["col10"], c10)))
    assert all(tuple(compare(a, b) for a, b in zip(new_table["col11"], c11)))

    table.as_msgpack_file(tmpdir.join("reflections.mpack").strpath)
    new_table = flex.reflection_table.from_msgpack_file(
        tmpdir.join("reflections.mpack").strpath
    )
    assert new_table.is_consistent()
    assert new_table.nrows() == 10
    assert new_table.ncols() == 11
    assert all(tuple(a == b for a, b in zip(new_table["col1"], c1)))
    assert all(tuple(a == b for a, b in zip(new_table["col2"], c2)))
    assert all(tuple(a == b for a, b in zip(new_table["col3"], c3)))
    assert all(tuple(a == b for a, b in zip(new_table["col4"], c4)))
    assert all(tuple(a == b for a, b in zip(new_table["col5"], c5)))
    assert all(tuple(a == b for a, b in zip(new_table["col6"], c6)))
    assert all(tuple(a == b for a, b in zip(new_table["col7"], c7)))
    assert all(tuple(a == b for a, b in zip(new_table["col8"], c8)))
    assert all(tuple(a == b for a, b in zip(new_table["col9"], c9)))
    assert all(tuple(a == b for a, b in zip(new_table["col10"], c10)))
    assert all(tuple(compare(a, b) for a, b in zip(new_table["col11"], c11)))


def test_experiment_identifiers():
    from dxtbx.model import ExperimentList, Experiment

    table = flex.reflection_table()
    table["id"] = flex.int([0, 1, 2, 3])

    table.assert_experiment_identifiers_are_consistent()

    identifiers = table.experiment_identifiers()
    identifiers[0] = "abcd"
    identifiers[1] = "efgh"
    identifiers[2] = "ijkl"
    identifiers[3] = "mnop"

    assert identifiers[0] == "abcd"
    assert identifiers[1] == "efgh"
    assert identifiers[2] == "ijkl"
    assert identifiers[3] == "mnop"

    for k, v in identifiers:
        if k == 0:
            assert v == "abcd"
        if k == 1:
            assert v == "efgh"
        if k == 2:
            assert v == "ijkl"
        if k == 3:
            assert v == "mnop"

    assert tuple(identifiers.keys()) == (0, 1, 2, 3)
    assert tuple(identifiers.values()) == ("abcd", "efgh", "ijkl", "mnop")

    table.assert_experiment_identifiers_are_consistent()

    experiments = ExperimentList()
    experiments.append(Experiment(identifier="abcd"))
    experiments.append(Experiment(identifier="efgh"))
    experiments.append(Experiment(identifier="ijkl"))
    experiments.append(Experiment(identifier="mnop"))

    table.assert_experiment_identifiers_are_consistent()

    experiments = ExperimentList()
    experiments.append(Experiment(identifier="abcd"))
    experiments.append(Experiment(identifier="efgh"))
    experiments.append(Experiment(identifier="ijkl"))
    experiments.append(Experiment(identifier="mnop"))
    experiments[3].identifier = "ijkl"

    with pytest.raises(AssertionError):
        table.assert_experiment_identifiers_are_consistent(experiments)

    experiments[2].identifier = "mnop"
    table.assert_experiment_identifiers_are_consistent(experiments)

    identifiers = table.experiment_identifiers()
    identifiers[0] = "abcd"
    identifiers[1] = "efgh"
    identifiers[2] = "ijkl"
    identifiers[3] = "ijkl"

    with pytest.raises(AssertionError):
        table.assert_experiment_identifiers_are_consistent()

    identifiers[3] = "mnop"

    import six.moves.cPickle as pickle

    pickled = pickle.dumps(table)
    table2 = pickle.loads(pickled)

    id1 = table.experiment_identifiers()
    id2 = table2.experiment_identifiers()

    for i in id1.keys():
        assert id1[i] == id2[i]

    other_table = flex.reflection_table()
    other_table["id"] = flex.int([3, 4])

    table.assert_experiment_identifiers_are_consistent()

    packed = table.as_msgpack()
    table2 = table.from_msgpack(packed)

    id1 = table.experiment_identifiers()
    id2 = table2.experiment_identifiers()

    for i in id1.keys():
        assert id1[i] == id2[i]

    other_table = flex.reflection_table()
    other_table["id"] = flex.int([3, 4])

    table.assert_experiment_identifiers_are_consistent()

    identifiers = other_table.experiment_identifiers()
    identifiers[3] = "mnop"
    identifiers[4] = "qrst"

    table.extend(other_table)

    assert len(table.experiment_identifiers()) == 5
    assert table.experiment_identifiers()[0] == "abcd"
    assert table.experiment_identifiers()[1] == "efgh"
    assert table.experiment_identifiers()[2] == "ijkl"
    assert table.experiment_identifiers()[3] == "mnop"
    assert table.experiment_identifiers()[4] == "qrst"

    assert len(table.experiment_identifiers()) == 5
    assert table.experiment_identifiers()[0] == "abcd"
    assert table.experiment_identifiers()[1] == "efgh"
    assert table.experiment_identifiers()[2] == "ijkl"
    assert table.experiment_identifiers()[3] == "mnop"
    assert table.experiment_identifiers()[4] == "qrst"


def test_select_remove_on_experiment_identifiers():

    table = flex.reflection_table()
    table["id"] = flex.int([0, 1, 2, 3])

    experiments = ExperimentList()
    experiments.append(Experiment(identifier="abcd"))
    experiments.append(Experiment(identifier="efgh"))
    experiments.append(Experiment(identifier="ijkl"))
    experiments.append(Experiment(identifier="mnop"))
    table.experiment_identifiers()[0] = "abcd"
    table.experiment_identifiers()[1] = "efgh"
    table.experiment_identifiers()[2] = "ijkl"
    table.experiment_identifiers()[3] = "mnop"

    table.assert_experiment_identifiers_are_consistent(experiments)

    table = table.remove_on_experiment_identifiers(["efgh"])
    del experiments[1]
    table.assert_experiment_identifiers_are_consistent(experiments)

    assert list(table.experiment_identifiers().keys()) == [0, 2, 3]
    assert list(table.experiment_identifiers().values()) == ["abcd", "ijkl", "mnop"]

    table = table.select_on_experiment_identifiers(["abcd", "mnop"])
    del experiments[1]  # now ijkl
    table.assert_experiment_identifiers_are_consistent(experiments)
    assert list(table.experiment_identifiers().keys()) == [0, 3]
    assert list(table.experiment_identifiers().values()) == ["abcd", "mnop"]

    # reset 'id' column such that they are numbered 0 .. n-1
    table.reset_ids()
    table.assert_experiment_identifiers_are_consistent(experiments)
    assert list(table.experiment_identifiers().keys()) == [0, 1]
    assert list(table.experiment_identifiers().values()) == ["abcd", "mnop"]
    # test that the function doesn't fail if no identifiers set
    table1 = copy.deepcopy(table)
    for k in table1.experiment_identifiers().keys():
        del table1.experiment_identifiers()[k]
    table1.reset_ids()
    assert list(table1.experiment_identifiers().keys()) == []

    # Test exception is raised if bad choice
    with pytest.raises(KeyError):
        table.remove_on_experiment_identifiers(["efgh"])
    with pytest.raises(KeyError):
        table.select_on_experiment_identifiers(["efgh"])

    table = flex.reflection_table()
    table["id"] = flex.int([0, 1, 2, 3])
    # Test exception is raised if identifiers map not set
    with pytest.raises(KeyError):
        table.remove_on_experiment_identifiers(["efgh"])
    with pytest.raises(KeyError):
        table.select_on_experiment_identifiers(["abcd", "mnop"])


def test_as_miller_array():
    table = flex.reflection_table()
    table["intensity.1.value"] = flex.double([1.0, 2.0, 3.0])
    table["intensity.1.variance"] = flex.double([0.25, 1.0, 4.0])
    table["miller_index"] = flex.miller_index([(1, 0, 0), (2, 0, 0), (3, 0, 0)])

    crystal = Crystal(
        real_space_a=(10, 0, 0),
        real_space_b=(0, 11, 0),
        real_space_c=(0, 0, 12),
        space_group=sgtbx.space_group_info("P 222").group(),
    )
    experiment = Experiment(crystal=crystal)

    iobs = table.as_miller_array(experiment, intensity="1")
    assert list(iobs.data()) == list(table["intensity.1.value"])
    assert list(iobs.sigmas()) == list(table["intensity.1.variance"] ** 0.5)

    with pytest.raises(Sorry):
        _ = table.as_miller_array(experiment, intensity="2")
    table["intensity.2.value"] = flex.double([1.0, 2.0, 3.0])
    with pytest.raises(Sorry):
        _ = table.as_miller_array(experiment, intensity="2")


def test_map_centroids_to_reciprocal_space(dials_regression):
    data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")
    pickle_path = os.path.join(data_dir, "full.pickle")
    sweep_path = os.path.join(data_dir, "datablock_orig.json")

    refl = flex.reflection_table.from_pickle(pickle_path)
    datablock = load.datablock(sweep_path, check_format=False)[0]
    imageset = datablock.extract_imagesets()[0]
    detector = imageset.get_detector()
    scan = imageset.get_scan()
    beam = imageset.get_beam()
    goniometer = imageset.get_goniometer()

    # check mm values not in
    assert "xyzobs.mm.value" not in refl

    refl.centroid_px_to_mm(detector, scan=scan)

    for k in ("xyzobs.mm.value", "xyzobs.mm.variance"):
        assert k in refl

    assert refl["xyzobs.mm.value"][0] == pytest.approx(
        (199.43400000000003, 11.908133333333334, 1.4324789835743459)
    )
    assert refl["xyzobs.mm.variance"][0] == pytest.approx(
        (0.0035346345381526106, 0.0029881028112449803, 5.711576621000785e-07)
    )

    refl.map_centroids_to_reciprocal_space(detector, beam, goniometer=goniometer)

    for k in ("s1", "rlp"):
        assert k in refl

    assert refl["s1"][0] == pytest.approx(
        (-0.035321308540942425, 0.6030297672949761, -0.8272574664632307)
    )
    assert refl["rlp"][0] == pytest.approx(
        (-0.035321308540942425, 0.27833194706770875, -0.5700990597173606)
    )

    # select only those centroids on the first image
    sel = refl["xyzobs.px.value"].parts()[2] < 1
    refl1 = refl.select(sel)
    del refl1["xyzobs.mm.value"], refl1["xyzobs.mm.variance"], refl1["s1"], refl1["rlp"]

    # pretend this is a still and hence no scan or goniometer
    refl1.centroid_px_to_mm(detector, scan=None)
    refl1.map_centroids_to_reciprocal_space(detector, beam, goniometer=None)

    assert refl1["s1"][0] == pytest.approx(
        (-0.035321308540942425, 0.6030297672949761, -0.8272574664632307)
    )
    # numbers for rlp are different to above since for the goniometer case the
    # starting angle of the first image is non-zero, so the rlps are rotated back
    # to zero degrees
    assert refl1["rlp"][0] == pytest.approx(
        (-0.035321308540942425, 0.6030297672949761, 0.19707031842793443)
    )


def test_calculate_entering_flags(dials_regression):
    data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")
    pickle_path = os.path.join(data_dir, "full.pickle")
    experiments_path = os.path.join(data_dir, "experiments_import.json")

    refl = flex.reflection_table.from_pickle(pickle_path)
    experiments = load.experiment_list(experiments_path, check_format=False)
    experiment = experiments[0]
    detector = experiment.detector
    scan = experiment.scan
    beam = experiment.beam
    goniometer = experiment.goniometer

    refl.centroid_px_to_mm(detector, scan=scan)
    refl.map_centroids_to_reciprocal_space(detector, beam, goniometer=goniometer)
    refl.calculate_entering_flags(experiments)
    assert "entering" in refl
    flags = refl["entering"]
    assert flags.count(True) == 58283
    assert flags.count(False) == 57799
