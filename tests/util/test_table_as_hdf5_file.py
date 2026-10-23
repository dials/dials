from __future__ import annotations

import os

import h5py
import hdf5plugin

from dials.algorithms.shoebox import MaskCode
from dials.array_family import flex
from dials.model.data import Shoebox
from dials.util.table_as_hdf5_file import HDF5TableFile


def test_table_as_hdf5_file_no_sbox(dials_data, tmp_path):
    data = dials_data("l_cysteine_4_sweeps_scaled") / "scaled_20_25.refl"
    table = flex.reflection_table.from_file(data)

    os.environ["DIALS_USE_H5"] = "1"
    table.as_file(tmp_path / "scaled.refl")

    data = h5py.File(tmp_path / "scaled.refl", "r")
    # try reading the data to check it was written as h5 as a single group
    assert len(data["dials"]["processing"].items()) == 1

    h5_refls = flex.reflection_table.from_file(tmp_path / "scaled.refl")

    def simple_test_equal(t1, t2):
        assert t1["intensity.sum.value"] == t2["intensity.sum.value"]
        assert list(t1.keys()) == list(t2.keys())

    simple_test_equal(h5_refls, table)

    # now test splitting and saving manually, both together
    split = table.split_by_experiment_id()
    with HDF5TableFile(tmp_path / "scaled2.refl", "w") as handle:
        handle.add_tables(split)

    data = h5py.File(tmp_path / "scaled2.refl", "r")
    # try reading the two experiments to check it was written as h5
    dset0 = data["dials"]["processing"]["group_0"]
    dset1 = data["dials"]["processing"]["group_1"]
    assert dset0["intensity.sum.value"].shape == (4417,)
    assert dset1["intensity.sum.value"].shape == (5555,)

    # now test saving manually, one at a time
    with HDF5TableFile(tmp_path / "scaled3.refl", "w") as handle:
        handle.add_tables([split[0]])
        handle.add_tables([split[1]])

    data = h5py.File(tmp_path / "scaled3.refl", "r")
    # try reading the two experiments to check it was written as h5
    dset0 = data["dials"]["processing"]["group_0"]
    dset1 = data["dials"]["processing"]["group_1"]
    assert dset0["intensity.sum.value"].shape == (4417,)
    assert dset1["intensity.sum.value"].shape == (5555,)

    with HDF5TableFile(tmp_path / "scaled3.refl", "r") as handle:
        tables = handle.get_tables()

    simple_test_equal(tables[0], split[0])
    simple_test_equal(tables[1], split[1])


def test_table_as_hdf5_file_compression(tmp_path):
    shoeboxes = flex.shoebox()
    for n in range(5):
        sbox = Shoebox(0, (0, 4, 0, 3, n, n + 2))
        sbox.allocate_data()
        sbox.allocate_background()
        for i in range(sbox.data.size()):
            sbox.data[i] = float(i % 7 + n)
            sbox.background[i] = 1.5
            sbox.mask[i] = MaskCode.Valid | (MaskCode.Foreground if i % 3 else 0)
        shoeboxes.append(sbox)
    table = flex.reflection_table()
    table["shoebox"] = shoeboxes
    table["id"] = flex.int(len(table), 0)
    table["bbox"] = shoeboxes.bounding_boxes()
    table["xyzobs.px.value"] = flex.vec3_double(len(table), (1.5, 2.5, 3.5))
    table["flags"] = flex.size_t(len(table), 32)
    table["entering"] = flex.bool(len(table), True)
    table["label"] = flex.std_string(len(table), "a")

    table.as_hdf5(tmp_path / "compressed.refl")

    # Numeric columns and the shoebox arrays use Bitshuffle with LZ4
    with h5py.File(tmp_path / "compressed.refl", "r") as handle:
        group = handle["dials"]["processing"]["group_0"]
        datasets = [group[k] for k in table.keys() if k not in ("shoebox", "label")]
        datasets.extend(group["shoebox"].values())
        assert len(datasets) == 10
        for dataset in datasets:
            assert str(hdf5plugin.BSHUF_ID) in dataset._filters, dataset.name

    # and round trip exactly
    new_table = flex.reflection_table.from_file(tmp_path / "compressed.refl")
    assert list(new_table.keys()) == list(table.keys())
    for new_array, array in zip(
        new_table["shoebox"].get_shoebox_data_arrays(),
        table["shoebox"].get_shoebox_data_arrays(),
    ):
        assert new_array.all_eq(array)
    assert list(new_table["shoebox"].bounding_boxes()) == list(
        table["shoebox"].bounding_boxes()
    )
    for key in ("id", "bbox", "xyzobs.px.value", "flags", "entering", "label"):
        assert list(new_table[key]) == list(table[key]), key
