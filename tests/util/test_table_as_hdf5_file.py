from __future__ import annotations

import os

import h5py

from dials.array_family import flex
from dials.util.table_as_hdf5_file import HDF5TableFile


def test_table_as_hdf5_file_no_sbox(dials_data, tmp_path):
    data = dials_data("l_cysteine_4_sweeps_scaled", pathlib=True) / "scaled_20_25.refl"
    table = flex.reflection_table.from_file(data)

    os.environ["DIALS_USE_H5"] = "1"
    table.as_file(tmp_path / "scaled.refl")

    data = h5py.File(tmp_path / "scaled.refl", "r")
    # try reading the data to check it was written as h5
    assert len(data["entry"]["data_processing"].items()) == 1

    h5_refls = flex.reflection_table.from_file(tmp_path / "scaled.refl")

    def simple_test_equal(t1, t2):
        assert t1["intensity.sum.value"] == t2["intensity.sum.value"]
        assert list(t1.keys()) == list(t2.keys())

    simple_test_equal(h5_refls, table)

    # now test splitting and saving manually, both together
    split = table.split_by_experiment_id()
    with HDF5TableFile(tmp_path / "scaled2.refl", "w") as handle:
        handle.set_reflections(split)

    data = h5py.File(tmp_path / "scaled2.refl", "r")
    # try reading the two experiments to check it was written as h5
    dset0 = data["entry"]["data_processing"]["0"]
    dset1 = data["entry"]["data_processing"]["1"]
    assert dset0.attrs["num_reflections"] == 4417
    assert dset1.attrs["num_reflections"] == 5555

    # now test saving manually, one at a time
    with HDF5TableFile(tmp_path / "scaled3.refl", "w") as handle:
        handle.set_reflections([split[0]])
        handle.set_reflections([split[1]])

    data = h5py.File(tmp_path / "scaled3.refl", "r")
    # try reading the two experiments to check it was written as h5
    dset0 = data["entry"]["data_processing"]["0"]
    dset1 = data["entry"]["data_processing"]["1"]
    assert dset0.attrs["num_reflections"] == 4417
    assert dset1.attrs["num_reflections"] == 5555

    with HDF5TableFile(tmp_path / "scaled3.refl", "r") as handle:
        tables = handle.get_reflections()

    simple_test_equal(tables[0], split[0])
    simple_test_equal(tables[1], split[1])
