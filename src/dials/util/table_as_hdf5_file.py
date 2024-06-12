from __future__ import annotations

from typing import Dict, List, Optional

import h5py
import hdf5plugin
import numpy as np

from dxtbx import flumpy

import dials_array_family_flex_ext
from dials.array_family import flex


class ReflectionListEncoder(object):
    """Encoder for the reflection data."""

    @staticmethod
    def encode(
        reflections: List[flex.reflection_table],
        handle: h5py.File,
    ) -> None:
        """Encode the list of reflection tables into a per-experiment hdf5 group."""

        # Create the reflection data group
        group = handle.create_group("entry/data_processing", track_order=True)

        for table in reflections:
            identifier_map = dict(table.experiment_identifiers())
            assert len(identifier_map) == 1
            identifier = list(identifier_map.values())[0]
            this_group = group.create_group(identifier)
            this_group.attrs["num_reflections"] = table.size()
            ReflectionListEncoder.encode_columns(
                this_group,
                table,
                ignore=["id"],
            )

    @staticmethod
    def encode_columns(
        group: h5py.Group,
        table: flex.reflection_table,
        ignore: Optional[List[str]] = None,
    ) -> None:
        """Encode the columns of a reflection table."""

        for key, data in table.cols():
            if ignore and key in ignore:
                continue
            if isinstance(data, flex.shoebox):
                ReflectionListEncoder.encode_shoebox(group, data, key)
            elif isinstance(data, flex.int6):
                s = data.size()
                this_data = flumpy.to_numpy(data.as_int()).reshape(s, 6)
                group.create_dataset(
                    key, data=this_data, shape=this_data.shape, dtype=this_data.dtype
                )
            elif isinstance(data, flex.std_string):
                this_data = [s.encode("utf-8") for s in data]
                group.create_dataset(key, data=this_data, shape=(len(this_data),))
            else:
                this_data = flumpy.to_numpy(data)
                group.create_dataset(
                    key, data=this_data, shape=this_data.shape, dtype=this_data.dtype
                )

    @staticmethod
    def encode_shoebox(group: h5py.Group, data: flex.shoebox, key: str):
        """Encode a column of shoeboxes."""
        sbdata, bg, mask = data.get_shoebox_data_arrays()
        panel = data.panels()
        bbox = data.bounding_boxes()
        data = flumpy.to_numpy(sbdata)
        lz4 = hdf5plugin.LZ4()
        sbox_group = group.create_group(key)
        sbox_group.create_dataset(
            "shoebox_data",
            data=data,
            shape=data.shape,
            dtype=data.dtype,
            compression=lz4,
        )
        data = flumpy.to_numpy(bg)
        sbox_group.create_dataset(
            "shoebox_background",
            data=data,
            shape=data.shape,
            dtype=data.dtype,
            compression=lz4,
        )
        data = flumpy.to_numpy(mask)
        sbox_group.create_dataset(
            "shoebox_mask",
            data=data,
            shape=data.shape,
            dtype=data.dtype,
            compression=lz4,
        )
        s = bbox.size()
        this_data = flumpy.to_numpy(bbox.as_int()).reshape(s, 6)
        sbox_group.create_dataset(
            "bbox",
            data=this_data,
            shape=this_data.shape,
            dtype=this_data.dtype,
            compression=lz4,
        )
        data = flumpy.to_numpy(panel)
        sbox_group.create_dataset(
            "panel", data=data, shape=data.shape, dtype=data.dtype, compression=lz4
        )


class ReflectionListDecoder(object):
    """Decoder for the reflection data."""

    @staticmethod
    def decode(handle: h5py.File) -> List[flex.reflection_table]:
        """Decode the data to a list of reflection tables."""

        # Get the group containing the reflection data
        g = handle["entry/data_processing"]

        # Create the list of reflection tables
        tables = []
        for i, (name, dataset) in enumerate(g.items()):
            if "num_reflections" not in dataset.attrs:
                raise RuntimeError(
                    "Unable to understand file as h5 reflection data (no num_reflections attribute)"
                )
            table = flex.reflection_table(int(dataset.attrs["num_reflections"]))
            table.experiment_identifiers()[i] = str(name)
            table["id"] = flex.int(table.size(), i)

            for key in dataset:
                if isinstance(dataset[key], h5py.Group):
                    # Decode all the shoebox data
                    shoebox_arrays: Dict[str, flumpy.FlexArray] = {}
                    names = [
                        "shoebox_data",
                        "shoebox_background",
                        "shoebox_mask",
                        "panel",
                        "bbox",
                    ]
                    for k in dataset[key].keys():
                        if k in names[:-1]:
                            shoebox_arrays[k] = flumpy.from_numpy(
                                np.array(dataset[key][k])
                            )
                        elif k == "bbox":
                            shoebox_arrays[k] = flex.int6(
                                flumpy.from_numpy(np.array(dataset[key][k]).flatten())
                            )
                        else:
                            raise RuntimeError(
                                f"Unrecognised elements {k} in {dataset[key]}"
                            )
                    if all(n in shoebox_arrays.keys() for n in names):
                        table[key] = flex.shoebox(
                            shoebox_arrays["panel"],
                            shoebox_arrays["bbox"],
                            allocate=True,
                        )
                        dials_array_family_flex_ext.ShoeboxExtractFromData(
                            table[key],
                            shoebox_arrays["shoebox_data"],
                            shoebox_arrays["shoebox_background"],
                            shoebox_arrays["shoebox_mask"],
                        )
                else:
                    table[key] = ReflectionListDecoder.convert_array(dataset[key])

            tables.append(table)

        # Return the list of reflection tables
        return tables

    @staticmethod
    def convert_array(data: np.array) -> flumpy.FlexArray:
        # Must allow that the data were written by a program outside of DIALS, so no special
        # knowledge of flex types should be assumed.
        if len(data.shape) == 2:
            if data.shape[1] == 3 and np.issubdtype(data.dtype, np.integer):
                # there is no such thing as flex.vec3_int, so this must be a miller index
                new = flumpy.miller_index_from_numpy(np.array(data))
            elif data.shape[1] == 3 or data.shape[1] == 2:  # vec3 or vec2 double
                new = flumpy.vec_from_numpy(np.array(data))
            elif data.shape[1] == 6 and np.issubdtype(data.dtype, np.integer):
                new = flex.int6(flumpy.from_numpy(np.array(data).flatten()))
                # N.B. flumpy could support this - but int6 only currently defined in dials, so would have to
                # move that to dxtbx first.
            else:
                raise RuntimeError(
                    f"Unrecognised 2D data dimensions/types: {data.shape} (expected shape (N,m)) where m is 2,3 or 6"
                )
        elif len(data.shape) == 1:
            if data.dtype == np.dtype("O"):  # "object type", for flex.std_string
                new = flex.std_string([s.decode("utf-8") for s in data])
            else:
                new = flumpy.from_numpy(np.array(data))
        elif len(data.shape) == 3:
            if data.shape[1] == 3 and data.shape[2] == 3:
                # it's a mat3 double
                new = flumpy.mat3_from_numpy(np.array(data))
            else:
                raise RuntimeError(
                    f"Unrecognised 3D data dimensions: {data.shape} (expected shape (N,3,3))"
                )

        else:
            raise RuntimeError(
                f"Unrecognised data dimensions: {data.shape}. Only 1D, 2D or 3D data supported."
            )
        return new


class HDF5TableFile:
    """
    Interface to on-disk representation of reflection data in hdf5 format.

    Note that for a multi-experiment tables, data is saved in a per-experiment
    table.
    """

    def __init__(self, filename: str, mode="w") -> None:
        """Open the file with the given mode."""
        self._handle = h5py.File(filename, mode)

    def close(self) -> None:
        """Close the file."""
        self._handle.close()
        del self._handle

    def set_data(
        self, data: List[flex.reflection_table], encoder: ReflectionListEncoder
    ) -> None:
        """Set the model data using the supplied encoder."""
        encoder.encode(data, self._handle)

    def get_data(self, decoder: ReflectionListDecoder) -> List[flex.reflection_table]:
        """Get the model data using the supplied decoder."""
        return decoder.decode(self._handle)

    def set_reflections(self, reflections: List[flex.reflection_table]) -> None:
        """Set the reflection data."""
        self.set_data(reflections, ReflectionListEncoder())

    def get_reflections(self) -> List[flex.reflection_table]:
        """Get the reflection data."""
        return self.get_data(ReflectionListDecoder())
