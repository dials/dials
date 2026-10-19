"""
Writing a reflection table with many experiment identifiers to the HDF5
.refl format fails with

    OSError: Unable to synchronously create attribute
             (object header message is too large)

dials.util.table_as_hdf5_file stores the experiment identifiers as an HDF5
attribute ("identifiers", a list of variable-length strings) on the data
group. With the default HDF5 file format, attributes live in the object
header, whose messages are limited to 64 KiB. Each variable-length string
costs a 16-byte heap reference there, so the write fails at exactly 4092
identifiers regardless of string length (4091 works). The msgpack writer has
no such limit.

Hit in practice by dials.combine_experiments on ~4500 stills experiments.

Possible fixes (both verified with 10000 identifiers, h5py 3.16/HDF5 1.14.6):
open the file with h5py.File(..., libver="latest") so attributes use dense
storage, or store the identifiers as a dataset instead of an attribute.
"""

from __future__ import annotations

import uuid

import pytest

from dials.array_family import flex

# Largest number of identifiers the attribute can hold with the default
# HDF5 file format; one more overflows the 64 KiB object header message.
MAX_IDENTIFIERS_IN_ATTRIBUTE = 4091


def _table_with_identifiers(n_expts, refl_per_expt=2):
    n = n_expts * refl_per_expt
    table = flex.reflection_table()
    table["id"] = flex.int([i // refl_per_expt for i in range(n)])
    table["panel"] = flex.size_t(n, 0)
    table["xyzobs.px.value"] = flex.vec3_double(n, (1.0, 2.0, 0.0))
    identifiers = table.experiment_identifiers()
    for i in range(n_expts):
        identifiers[i] = str(uuid.uuid4())
    return table


def _check_roundtrip(table, filename, writer):
    writer(table, filename)
    back = flex.reflection_table.from_file(filename)
    assert len(back) == len(table)
    assert dict(back.experiment_identifiers()) == dict(table.experiment_identifiers())


def test_hdf5_identifiers_at_attribute_limit(tmp_path):
    table = _table_with_identifiers(MAX_IDENTIFIERS_IN_ATTRIBUTE)
    _check_roundtrip(table, tmp_path / "ok.refl", lambda t, fn: t.as_hdf5(fn))


@pytest.mark.xfail(
    raises=OSError,
    strict=True,
    reason="identifiers attribute exceeds the 64 KiB HDF5 object header limit",
)
def test_hdf5_identifiers_above_attribute_limit(tmp_path):
    table = _table_with_identifiers(MAX_IDENTIFIERS_IN_ATTRIBUTE + 1)
    _check_roundtrip(table, tmp_path / "fail.refl", lambda t, fn: t.as_hdf5(fn))


def test_msgpack_identifiers_above_attribute_limit(tmp_path):
    table = _table_with_identifiers(MAX_IDENTIFIERS_IN_ATTRIBUTE + 1)
    _check_roundtrip(table, tmp_path / "ok.refl", lambda t, fn: t.as_msgpack_file(fn))
