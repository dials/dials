from __future__ import annotations

from dials.array_family import flex
from dials.model.data import Shoebox


def gen_shoebox():
    shoebox = Shoebox(0, (0, 4, 0, 3, 0, 1))
    shoebox.allocate()
    for k in range(1):
        for j in range(3):
            for i in range(4):
                shoebox.data[k, j, i] = i + 2 * j + 3 * k + 0.1
                shoebox.mask[k, j, i] = i % 2
                shoebox.background[k, j, i] = i * j + 0.2
    return shoebox


def table_and_columns():
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
    c11 = [gen_shoebox() for _ in range(10)]

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
    return table, [c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11]


if __name__ == "__main__":
    table, columns = table_and_columns()
    c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11 = columns
    table["id"] = flex.int(table.size(), 0)
    table.experiment_identifiers()[0] = "test"
    ds = table.to_xarray()
    print(ds)
