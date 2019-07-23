#!/usr/bin/env python
#
# dials.util.image.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division, print_function


class reader:
    """A class to read the CBF files used in DIALS"""

    def read_file(self, filename):
        """Read the CBF file"""
        import pycbf

        self.cbf_handle = pycbf.cbf_handle_struct()
        self.cbf_handle.read_file(filename, pycbf.MSG_DIGEST)
        self.cbf_handle.rewind_datablock()

    def get_data(self):
        """Get the gain array from the file"""
        import numpy
        from scitbx.array_family import flex

        # Select the first datablock and rewind all the categories
        self.cbf_handle.select_datablock(0)
        self.cbf_handle.select_category(0)
        self.cbf_handle.select_column(2)
        self.cbf_handle.select_row(0)

        # Check the type of the element to ensure it's a binary
        # otherwise raise an exception
        if "bnry" in self.cbf_handle.get_typeofvalue():
            # Read the image data into an array
            image_string = self.cbf_handle.get_integerarray_as_string()
            image = flex.int(numpy.fromstring(image_string, numpy.int32))

            # Get the array parameters
            parameters = self.cbf_handle.get_integerarrayparameters_wdims()
            image_size = (parameters[10], parameters[9])

            # Resize the image
            image.reshape(flex.grid(*image_size))
        else:
            raise TypeError("Can't find image")

        # Return the image
        return image


if __name__ == "__main__":
    import sys

    handle = reader()
    handle.read_file(sys.argv[1])
    image = handle.get_data()
