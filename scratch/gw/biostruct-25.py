from __future__ import division
# copied from James's code...

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
  import numpy

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

def open_file_return_array(filename):
  import pycbf
  import numpy

  # open file
  cbf_handle = pycbf.cbf_handle_struct()
  cbf_handle.read_file(filename, pycbf.MSG_DIGEST)
  cbf_handle.rewind_datablock()

  # find right datablock
  cbf_handle.select_datablock(0)
  cbf_handle.select_category(0)
  cbf_handle.select_column(2)
  cbf_handle.select_row(0)

  type = cbf_handle.get_typeofvalue()
  assert (type.find('bnry') > -1)

  # read and reshape
  image = numpy.fromstring(cbf_handle.get_integerarray_as_string(),
                           numpy.int32)
  parameters = cbf_handle.get_integerarrayparameters_wdims()
  image.shape = (parameters[10], parameters[9])

  return image

def plot_image(image):
  from matplotlib import pyplot
  pyplot.imshow(image)
  pyplot.savefig('biostruct-25.png')

  return

def parse_xds_xparm_scan_info(xparm_file):
  values = map(float, open(xparm_file).read().split())

  assert(len(values) == 42)

  img_start, osc_start, osc_range = values[:3]

  return img_start, osc_start, osc_range

def read_integrate_hkl_apply_corrections(
        int_hkl, x_crns_file, y_crns_file, xparm_file, image_file):
  '''Read X corrections and Y corrections to arrays; for record in
  int_hkl read x, y positions and compute r.m.s. deviations between
  observed and calculated centroids with and without the corrections.
  N.B. will only compute offsets for reflections with I/sigma > 10.'''

  from scitbx import matrix
  import math

  x_corrections = open_file_return_array(x_crns_file)
  y_corrections = open_file_return_array(y_crns_file)

  sumxx_orig = 0.0
  n_obs = 0
  sumxx_corr = 0.0

  # FIXME find code to run prediction of reflection lists...

  img_start, osc_start, osc_range = parse_xds_xparm_scan_info(xparm_file)

  from rstbx.cftbx.coordinate_frame_converter import \
      coordinate_frame_converter

  cfc = coordinate_frame_converter(xparm_file)

  fast = cfc.get_c('detector_fast')
  slow = cfc.get_c('detector_slow')
  normal = fast.cross(slow)

  origin = cfc.get_c('detector_origin')
  distance = origin.dot(normal)

  A = cfc.get_c('real_space_a')
  B = cfc.get_c('real_space_b')
  C = cfc.get_c('real_space_c')

  cell = (A.length(), B.length(), C.length(), B.angle(C, deg = True),
          C.angle(A, deg = True), A.angle(B, deg = True))

  u, b = cfc.get_u_b()
  ub = u * b
  s0 = - cfc.get_c('sample_to_source') / cfc.get('wavelength')
  axis = cfc.get_c('rotation_axis')

  xcal = []
  xobs = []
  xcor = []
  ycal = []
  yobs = []
  ycor = []

  for record in open(int_hkl):
    if record.startswith('!'):
      continue
    values = map(float, record.split())
    i, sigi = values[3], values[4]

    if sigi < 0:
      continue

    if i / sigi < 10:
      continue

    # FIXME also read the GXPARM.XDS file and recompute the reflection
    # position from the Miller indices and model of experimental geometry
    # as xc, yc, zc... Ah, smarter: use the zc from the file to provide the
    # rotation, not a full prediction: only then need to worry about
    # computing the intersection.

    hkl = map(int, map(round, values[:3]))

    xc, yc, zc = values[5:8]
    angle = (zc - img_start + 1) * osc_range + osc_start
    d2r = math.pi / 180.0

    s = (ub * hkl).rotate(axis, angle * d2r) + s0
    s = s.normalize()
    r = s * distance / s.dot(normal) - origin

    xc = r.dot(fast) / 0.172
    yc = r.dot(slow) / 0.172

    xo, yo, zo = values[12:15]

    xc_orig = matrix.col((xc, yc))
    xo_orig = matrix.col((xo, yo))

    n_obs += 1
    sumxx_orig += (xo_orig - xc_orig).dot()

    # get the correcions for this pixel... N.B. 1 + ... lost due to C-style
    # rather than Fortran-style arrays

    ix4 = int(round((xc - 2) / 4))
    iy4 = int(round((yc - 2) / 4))

    # hmm.... Fortran multi-dimensional array ordering... identified by
    # bounds error other way around...

    dx = 0.1 * x_corrections[iy4, ix4]
    dy = 0.1 * y_corrections[iy4, ix4]

    xc_corr = matrix.col((xc + dx, yc + dy))

    sumxx_corr += (xo_orig - xc_corr).dot()

    # Put coords in array if they're in frame zero
    if (zo >= 0 and zo < 1):
      xcal.append(xc)
      xobs.append(xo)
      xcor.append(xc + dx)

      ycal.append(yc)
      yobs.append(yo)
      ycor.append(yc + dy)


  if (image_file):
    import pycbf
    cbf_handle = pycbf.cbf_handle_struct()
    cbf_handle.read_file(image_file, pycbf.MSG_DIGEST)
    image_array = get_image(cbf_handle)

    # Correct coords for matplotlib convention
    xcal = [x - 0.5 for x in xcal]
    xobs = [x - 0.5 for x in xobs]
    xcor = [x - 0.5 for x in xcor]

    ycal = [y - 0.5 for y in ycal]
    yobs = [y - 0.5 for y in yobs]
    ycor = [y - 0.5 for y in ycor]

    # Manually selected a spot
#      xp = 533
#      yp = 342
#      ix4 = int(round((xp - 2) / 4))
#      iy4 = int(round((yp - 2) / 4))
#      dx = 0.1 * x_corrections[iy4, ix4]
#      dy = 0.1 * y_corrections[iy4, ix4]
#      xpcor = [xp + dx]
#      ypcor = [yp + dy]


    # Show the image with points plotted on it
    from matplotlib import pylab, cm
    pylab.imshow(image_array, cmap=cm.Greys_r, interpolation='nearest',
                 origin='bottom')
    pylab.scatter(xcal, ycal, c='r', marker='x')
    pylab.scatter(xobs, yobs, c='b', marker='x')
    pylab.scatter(xcor, ycor, c='y', marker='x')
    #pylab.scatter(xpcor, ypcor, c='b', marker='8')
    pylab.show()

  return math.sqrt(sumxx_orig / n_obs), math.sqrt(sumxx_corr / n_obs)

def read_spot_xds_apply_corrections(
        spot_xds, x_crns_file, y_crns_file, xparm_file):
  '''Read X corrections and Y corrections to arrays; for record in
  int_hkl read x, y positions and compute r.m.s. deviations between
  observed and calculated centroids with and without the corrections.
  N.B. will only compute offsets for reflections with I/sigma > 10.'''

  from scitbx import matrix
  import math

  x_corrections = open_file_return_array(x_crns_file)
  y_corrections = open_file_return_array(y_crns_file)

  sumxx_orig = 0.0
  n_obs = 0
  sumxx_corr = 0.0

  # FIXME find code to run prediction of reflection lists...

  img_start, osc_start, osc_range = parse_xds_xparm_scan_info(xparm_file)

  from rstbx.cftbx.coordinate_frame_converter import \
      coordinate_frame_converter

  cfc = coordinate_frame_converter(xparm_file)

  fast = cfc.get_c('detector_fast')
  slow = cfc.get_c('detector_slow')
  normal = fast.cross(slow)

  origin = cfc.get_c('detector_origin')
  distance = origin.dot(normal)

  A = cfc.get_c('real_space_a')
  B = cfc.get_c('real_space_b')
  C = cfc.get_c('real_space_c')

  cell = (A.length(), B.length(), C.length(), B.angle(C, deg = True),
          C.angle(A, deg = True), A.angle(B, deg = True))

  u, b = cfc.get_u_b()
  ub = u * b
  s0 = - cfc.get_c('sample_to_source') / cfc.get('wavelength')
  axis = cfc.get_c('rotation_axis')

  xcal = []
  xobs = []
  xcor = []
  ycal = []
  yobs = []
  ycor = []

  for record in open(spot_xds):
    values = map(float, record.split())

    hkl = map(int, map(round, values[-3:]))

    if hkl[0] == 0 and hkl[1] == 0 and hkl[2] == 0:
      continue

    xc, yc, zc = values[:3]
    angle = (zc - img_start + 1) * osc_range + osc_start
    d2r = math.pi / 180.0

    s = (ub * hkl).rotate(axis, angle * d2r) + s0
    s = s.normalize()
    r = s * distance / s.dot(normal) - origin

    xc = r.dot(fast) / 0.172
    yc = r.dot(slow) / 0.172

    xo, yo, zo = values[:3]

    xc_orig = matrix.col((xc, yc))
    xo_orig = matrix.col((xo, yo))

    n_obs += 1
    sumxx_orig += (xo_orig - xc_orig).dot()

    # get the correcions for this pixel... N.B. 1 + ... lost due to C-style
    # rather than Fortran-style arrays

    ix4 = int(round((xc - 2) / 4))
    iy4 = int(round((yc - 2) / 4))

    # hmm.... Fortran multi-dimensional array ordering... identified by
    # bounds error other way around...

    dx = 0.1 * x_corrections[iy4, ix4]
    dy = 0.1 * y_corrections[iy4, ix4]

    xc_corr = matrix.col((xc + dx, yc + dy))

    sumxx_corr += (xo_orig - xc_corr).dot()

  print n_obs
  return math.sqrt(sumxx_orig / n_obs), math.sqrt(sumxx_corr / n_obs)

if __name__ == '__main__':
  import sys

  orig, corr = read_spot_xds_apply_corrections(
      sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

  print orig, corr
