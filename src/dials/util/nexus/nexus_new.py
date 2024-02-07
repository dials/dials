from __future__ import annotations

import h5py
import numpy as np

from dxtbx import flumpy


class ReflectionListEncoder(object):
    """Encoder for the reflection data."""

    def encode(self, reflections, handle):
        """Encode the reflection data."""

        # Create the reflection data group
        group = handle.create_group("entry/data_processing", track_order=True)
        # ok this works if we don't have -1 ids i.e. after #2567 just remove for now.
        if -1 in reflections["id"]:
            sel = reflections["id"] == -1
            reflections.del_selected(sel)
        tables = reflections.split_by_experiment_id()

        for table in tables:
            table.clean_experiment_identifiers_map()
            assert len(table.experiment_identifiers().values()) == 1
            identifier = list(table.experiment_identifiers().values())[0]
            this_group = group.create_group(identifier)
            this_group.attrs["num_reflections"] = table.size()
            del table["id"]
            self.encode_columns(this_group, table)

    def encode_columns(self, group, table):
        """Encode a column of data."""
        from dials.array_family import flex

        for key, data in table.cols():
            if isinstance(data, flex.shoebox):
                self.encode_shoebox(group, table)
            elif isinstance(data, flex.int6):
                s = data.size()
                this_data = flumpy.to_numpy(data.as_int()).reshape(s, 6)
                group.create_dataset(
                    key, data=this_data, shape=this_data.shape, dtype=this_data.dtype
                )
            else:
                this_data = flumpy.to_numpy(data)
                group.create_dataset(
                    key, data=this_data, shape=this_data.shape, dtype=this_data.dtype
                )

    def encode_shoebox(self, group, table):
        """Encode a column of shoeboxes."""
        sbdata, bg, mask = table.get_shoebox_data_arrays()
        data = flumpy.to_numpy(sbdata)
        group.create_dataset(
            "shoebox_data", data=data, shape=data.shape, dtype=data.dtype
        )
        data = flumpy.to_numpy(bg)
        group.create_dataset(
            "shoebox_background", data=data, shape=data.shape, dtype=data.dtype
        )
        data = flumpy.to_numpy(mask)
        group.create_dataset(
            "shoebox_mask", data=data, shape=data.shape, dtype=data.dtype
        )


class ReflectionListDecoder(object):
    """Decoder for the reflection data."""

    def decode(self, handle):
        """Decode the reflection data."""
        import dials_array_family_flex_ext
        from dials.array_family import flex

        # Get the group containing the reflection data
        g = handle["entry/data_processing"]

        # Create the list of reflections
        tables = []
        for i, (name, dataset) in enumerate(g.items()):
            table = flex.reflection_table(int(dataset.attrs["num_reflections"]))
            table.experiment_identifiers()[i] = str(name)
            table["id"] = flex.int(table.size(), i)

            # Decode all the columns
            shoebox_data = None
            shoebox_background = None
            shoebox_mask = None

            for key in dataset:
                item = dataset[key]
                if key == "shoebox_data":
                    shoebox_data = flumpy.from_numpy(np.array(dataset[key]))
                elif key == "shoebox_background":
                    shoebox_background = flumpy.from_numpy(np.array(dataset[key]))
                elif key == "shoebox_mask":
                    shoebox_mask = flumpy.from_numpy(np.array(dataset[key]))
                else:
                    val = self._convert(key, item)
                    if val:
                        table[str(key)] = val

            if "panel" in table and "bbox" in table:
                table["shoebox"] = flex.shoebox(
                    table["panel"], table["bbox"], allocate=True
                )
                if shoebox_mask and shoebox_background and shoebox_data:
                    dials_array_family_flex_ext.ShoeboxExtractFromData(
                        table, shoebox_data, shoebox_background, shoebox_mask
                    )
            tables.append(table)
        table = flex.reflection_table.concat(tables)
        # Return the list of reflections
        return table

    def _convert(self, key, data):
        # Must allow that the data were written by a program outside of DIALS, so no special
        # knowledge of flex types should be assumed.
        from dials.array_family import flex

        if key == "miller_index":  # special
            assert len(data.shape) == 2 and data.shape[1] == 3
            new = flumpy.miller_index_from_numpy(np.array(data))
        elif len(data.shape) == 2:
            if data.shape[1] == 3 or data.shape[1] == 2:  # vec3 or vec2 double
                new = flumpy.vec_from_numpy(np.array(data))
            elif data.shape[1] == 6:
                new = flex.int6(flumpy.from_numpy(np.array(data).flatten()))
                # N.B. flumpy could support this - but int6 only currently defined in dials, so would have to
                # move that to dxtbx first.
            else:
                raise RuntimeError(
                    "Unrecognised 2D data dimensions (expected shape (N,m)) where m is 2,3 or 6"
                )
        else:
            new = flumpy.from_numpy(np.array(data))
        return new


class NewNexusFile:
    """Interface to Nexus file."""

    def __init__(self, filename, mode="w"):
        """Open the file with the given mode."""
        self._handle = h5py.File(filename, mode)

    def close(self):
        """Close the file."""
        self._handle.close()
        del self._handle

    def set_data(self, data, encoder):
        """Set the model data using the supplied encoder."""
        encoder.encode(data, self._handle)

    def get_data(self, decoder):
        """Get the model data using the supplied decoder."""
        return decoder.decode(self._handle)

    def set_reflections(self, reflections):
        """Set the reflection data."""
        self.set_data(reflections, ReflectionListEncoder())

    def get_reflections(self):
        """Get the reflection data."""
        return self.get_data(ReflectionListDecoder())
