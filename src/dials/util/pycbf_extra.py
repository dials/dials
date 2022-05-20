"""This module defines some useful header functions for CBF files."""

from __future__ import annotations

import numpy as np
import pycbf

from scitbx import matrix


def compute_central_rotation_matrix(gonio):
    """Compute the central rotation matrix

    :param gonio: The goniometer struct
    :returns: The central rotation matrix
    """
    x = gonio.rotate_vector(0.5, 1, 0, 0)
    y = gonio.rotate_vector(0.5, 0, 1, 0)
    z = gonio.rotate_vector(0.5, 0, 0, 1)
    R = matrix.sqr(x + y + z).transpose()
    return R


def get_image(cbf_handle, category="array_data", column="data", row=0, element=0):
    """Read an image from a CBF file

    This function is a bit of a hack - I'm not sure what the general structure
    of a CBF file is like but for the data I have, it works. Reads an image
    from the location specified in the CBF file, otherwise raises an exception.

    :param cbf_handle: The handle to the CBF file
    :param category: Category in which the image is contained
    :param column: Column in which the image is contained
    :param row: Row in which image is contained
    :param element: Element in which image is contained
    :returns: An array of image data
    """
    # Find the given category, column and row
    cbf_handle.find_category(category)
    cbf_handle.find_column(column)
    cbf_handle.select_row(row)

    # Check the type of the element to ensure it's a binary
    # otherwise raise an exception
    if "bnry" in cbf_handle.get_typeofvalue():

        # Read the image data into an array
        image_string = cbf_handle.get_integerarray_as_string()
        image = np.fromstring(image_string, np.int32)

        # Get the size of the image data (either 2d or 3d)
        image_size = cbf_handle.get_image_size(element)

        # Resize the image
        image.shape = image_size

    else:
        raise TypeError(f"{category}:{column}:{row}:{element} is not an image")

    # Return the image
    return image


def get_image_volume(cbf_paths):
    """Load the image volume from the list of cbf_paths. The list of paths is
    assumed to be is order from 1->n.

    :param cbf_paths: The list of cbf files
    :param width The width (xsize) of the volume
    :param height The height (ysize) of the volume
    :returns: The 3D volume array
    """
    # Read the first image and get the size
    cbf_handle = pycbf.cbf_handle_struct()
    cbf_handle.read_file(cbf_paths[0], pycbf.MSG_DIGEST)
    image = get_image(cbf_handle)
    height, width = image.shape

    # Initialise the image volume
    num_slices = len(cbf_paths)
    volume = np.zeros(shape=(num_slices, height, width), dtype=np.int32)
    volume[0, :, :] = image

    # For each CBF file, read the image and put into the image volume
    for i, filename in enumerate(cbf_paths[1:]):
        cbf_handle = pycbf.cbf_handle_struct()
        cbf_handle.read_file(filename, pycbf.MSG_DIGEST)
        volume[i + 1, :, :] = get_image(cbf_handle)

    # Return the image volume
    return volume


def search_for_image_volume(search_path):
    """Load the CBF image volume

    Args:
        search_path: The CBF file search path

    Returns:
        The image volume
    """
    from glob import glob

    # Load the CBF image volume
    cbf_path = sorted(glob(search_path))
    return get_image_volume(cbf_path)
