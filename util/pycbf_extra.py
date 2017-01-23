"""This module defines some useful header functions for CBF files."""
from __future__ import absolute_import, division

from scitbx import matrix
import pycbf
import numpy

def print_info(cbf_path):
  """Print out a load of data held in the CBF file.

  This is by no means a full list of the data contained in the file, it's
  mainly for debugging and development purposes. The data that will be
  printed is the following:
  - The number of categories and the name of each category
  - The number of rows and columns and the name of each column
  - The type of each element in each row/column element.

  :param cbf_path: The path to the cbf file

  """
  # Read the CBF file
  cbf_handle = pycbf.cbf_handle_struct()
  cbf_handle.read_file(cbf_path, pycbf.MSG_DIGEST)
  cbf_handle.rewind_datablock()

  # Select the first datablock and rewind all the categories
  cbf_handle.select_datablock(0)
  cbf_handle.rewind_category()

  # Count the number of categories and loop through them
  num_categories = cbf_handle.count_categories()
  for i in range(num_categories):

    # Select the ith category and print its name
    cbf_handle.select_category(i)
    category_name = cbf_handle.category_name()
    print "Category:", i, category_name

    # Count the number of rows and columns in the category
    # and print them
    num_rows = cbf_handle.count_rows()
    num_cols = cbf_handle.count_columns()
    print "\tNum (rows, cols)", (num_rows, num_cols)

    # Rewind the columns and print the name of each
    cbf_handle.rewind_column()
    for i in range(num_cols):
      cbf_handle.select_column(i)
      column_name = cbf_handle.column_name()
      print '\tColumn:', i, column_name

    # Loop through all rows and columns and print the
    # type of the data stored in that table element
    for j in range(num_rows):
      cbf_handle.select_row(j)
      cbf_handle.rewind_column()
      print '\t\tRow:', j, cbf_handle.get_value()
      for i in range(num_cols):
        cbf_handle.select_column(i)
        type_of_value = cbf_handle.get_typeofvalue()
        if type_of_value.find('dblq') > -1:
          value = cbf_handle.get_value()
        elif type_of_value.find('text') > -1:
          value = cbf_handle.get_value()
          value = '\n\t\t\t'.join(value.split('\n'))
        elif type_of_value.find('word') > -1:
          value = cbf_handle.get_value()
        elif type_of_value.find('sglq') > -1:
          value = cbf_handle.get_value()
        else:
          value = '...'
        print "\t\tColumn", i, "Type:", type_of_value, value


def get_beam_direction(cbf_handle):
  """Find the beam direction (why is this not simpler in pycbf?)

  :param cbf_handle: The cbf file handle
  :returns: The beam vector

  """
  cbf_handle.find_category('axis')
  cbf_handle.find_column('equipment')
  cbf_handle.find_row('source')
  beam_direction = []
  for j in range(3):
    cbf_handle.find_column('vector[%d]' % (j + 1))
    beam_direction.append(cbf_handle.get_doublevalue())

  B = - matrix.col(beam_direction).normalize()
  return B


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


def get_image(cbf_handle, category='array_data', column='data', row=0,
              element=0):
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
  type = cbf_handle.get_typeofvalue()
  if type.find('bnry') > -1:

    # Read the image data into an array
    image_string = cbf_handle.get_integerarray_as_string()
    image = numpy.fromstring(image_string, numpy.int32)

    # Get the size of the image data (either 2d or 3d)
    image_size = cbf_handle.get_image_size(element)

    # Resize the image
    image.shape = (image_size)

  else:
    raise TypeError('{0}:{1}:{2}:{3} is not an image'.format(
                        category, column, row, element))

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
  volume = numpy.zeros(shape=(num_slices, height, width), dtype=numpy.int32)
  volume[0,:,:] = image

  # For each CBF file, read the image and put into the image volume
  for i, filename in enumerate(cbf_paths[1:]):
    cbf_handle = pycbf.cbf_handle_struct()
    cbf_handle.read_file(filename, pycbf.MSG_DIGEST)
    volume[i+1,:,:] = get_image(cbf_handle)

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
  cbf_path = glob(search_path)
  cbf_path.sort()
  return get_image_volume(cbf_path)
