#!/usr/bin/env python
#
# nexus.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

from __future__ import absolute_import, division, print_function


class H5PYEncoder(object):
    """Encoder base class."""

    def encode(self, obj, handle):
        """Encode the object in the HDF5 file."""
        raise RuntimeError("Overload!")


class H5PYDecoder(object):
    """Decoder base class."""

    def decode(self, handle):
        """Decode the object from the HDF5 file."""
        raise RuntimeError("Overload!")


class ReflectionListEncoder(H5PYEncoder):
    """Encoder for the reflection data."""

    def encode(self, reflections, handle):
        """Encode the reflection data."""

        # Create the reflection data group
        group = handle.create_group("entry/data_processing")
        group.attrs["num_reflections"] = len(reflections)

        # Encode each column entry
        for key, data in reflections.cols():
            self.encode_column(group, key, data)

    def encode_column(self, group, key, data):
        """ Encode a column of data. """
        from dials.array_family import flex

        if isinstance(data, flex.shoebox):
            self.encode_shoebox(group, key, data)
        else:
            group[key] = list(data)
            group[key].attrs["flex_type"] = type(data).__name__

    def encode_shoebox(self, group, key, sb_data):
        """ Encode a column of shoeboxes. """
        shoebox = group.create_group("shoebox")
        shoebox.attrs["flex_type"] = type(sb_data).__name__
        data = shoebox.create_group("data")
        mask = shoebox.create_group("mask")
        background = shoebox.create_group("background")
        for i, sb in enumerate(sb_data):
            data["%d" % i] = sb.data.as_numpy_array()
            mask["%d" % i] = sb.mask.as_numpy_array()
            background["%d" % i] = sb.background.as_numpy_array()


class ReflectionListDecoder(H5PYDecoder):
    """Decoder for the reflection data."""

    def decode(self, handle):
        """Decode the reflection data."""
        from dials.array_family import flex

        # Get the group containing the reflection data
        g = handle["entry/data_processing"]

        # Create the list of reflections
        rl = flex.reflection_table(int(g.attrs["num_reflections"]))

        # Decode all the columns
        for key in g:
            item = g[key]
            name = item.attrs["flex_type"]
            if name == "shoebox":
                flex_type = getattr(flex, name)
                data = item["data"]
                mask = item["mask"]
                background = item["background"]
                col = flex_type(len(rl))
                for i in range(len(rl)):
                    col[i].data = flex.double(data["%d" % i].value)
                    col[i].mask = flex.int(mask["%d" % i].value)
                    col[i].background = flex.double(background["%d" % i].value)

            else:
                flex_type = getattr(flex, name)
                col = self.decode_column(flex_type, item)
            rl[str(key)] = col

        # Return the list of reflections
        return rl

    def decode_column(self, flex_type, data):
        """ Decode a column for various flex types. """
        from dials.array_family import flex

        if flex_type == flex.int:
            return flex_type(int(d) for d in list(data))
        elif flex_type == flex.double:
            return flex_type(float(d) for d in list(data))
        elif flex_type == flex.vec3_double:
            return flex_type(list(map(float, d)) for d in list(data))
        elif flex_type == flex.vec2_double:
            return flex_type(list(map(float, d)) for d in list(data))
        elif flex_type == flex.int6:
            return flex_type(list(map(int, d)) for d in list(data))
        elif flex_type == flex.miller_index:
            return flex_type(list(map(int, d)) for d in list(data))
        elif flex_type == flex.bool:
            return flex_type(bool(d) for d in list(data))
        elif flex_type == flex.size_t:
            return flex_type(int(d) for d in list(data))
        else:
            return flex_type(list(data))


class NexusFile(object):
    """Interface to Nexus file."""

    def __init__(self, filename, mode="a"):
        """Open the file with the given mode."""
        import h5py

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


if __name__ == "__main__":
    from dials.array_family import flex

    reflections = flex.reflection_table(
        [
            ("hkl", flex.miller_index(10)),
            ("s1", flex.vec3_double(10)),
            ("bbox", flex.int6(10)),
            ("id", flex.int(10)),
            ("shoebox", flex.shoebox(10)),
        ]
    )

    for i in range(10):
        reflections["shoebox"][i].data = flex.double(flex.grid(10, 10, 10))
        reflections["shoebox"][i].mask = flex.int(flex.grid(10, 10, 10))
        reflections["shoebox"][i].background = flex.double(flex.grid(10, 10, 10))

    for i in range(10):
        print(reflections["shoebox"][i].data.all())

    writer = NexusFile("temp.h5", "w")
    writer.set_reflections(reflections)
    writer.close()

    reader = NexusFile("temp.h5", "r")
    reflections = reader.get_reflections()
    reader.close()

    for r in reflections:
        print(r)
