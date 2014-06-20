def read_xds_calibration_file(calibration_file):
  '''Read XDS calibration file, return as flex array.'''

  from scitbx.array_family import flex
  from cbflib_adaptbx import uncompress, compress
  import binascii

  start_tag = binascii.unhexlify('0c1a04d5')

  data = open(calibration_file, 'rb').read()
  data_offset = data.find(start_tag) + 4
  cbf_header = data[:data_offset - 4]

  fast = 0
  slow = 0
  length = 0

  for record in cbf_header.split('\n'):
    if 'X-Binary-Size-Fastest-Dimension' in record:
      fast = int(record.split()[-1])
    elif 'X-Binary-Size-Second-Dimension' in record:
      slow = int(record.split()[-1])
    elif 'X-Binary-Number-of-Elements' in record:
      length = int(record.split()[-1])
    elif 'X-Binary-Size:' in record:
      size = int(record.split()[-1])

  assert(length == fast * slow)

  pixel_values = uncompress(packed = data[data_offset:data_offset + size],
                            fast = fast, slow = slow)

  return pixel_values

def run_job(executable, arguments = [], stdin = [], working_directory = None):
  '''Run a program with some command-line arguments and some input,
  then return the standard output when it is finished.'''

  import subprocess
  import os

  if working_directory is None:
    working_directory = os.getcwd()

  command_line = '%s' % executable
  for arg in arguments:
    command_line += ' "%s"' % arg

  popen = subprocess.Popen(command_line,
                           bufsize = 1,
                           stdin = subprocess.PIPE,
                           stdout = subprocess.PIPE,
                           stderr = subprocess.STDOUT,
                           cwd = working_directory,
                           universal_newlines = True,
                           shell = True,
                           env = os.environ)

  for record in stdin:
    popen.stdin.write('%s\n' % record)

  popen.stdin.close()

  output = []

  while True:
    record = popen.stdout.readline()
    if not record:
      break

    output.append(record)

  return output

def fit_plane(xyz_values):
  '''Derive equation of plane

  z = ax + by + c

  from list of values x, y, z; returns a, b, c.'''

  from scitbx import matrix

  x0 = sum([d[0] for d in xyz_values]) / len(xyz_values)
  y0 = sum([d[1] for d in xyz_values]) / len(xyz_values)
  z0 = sum([d[2] for d in xyz_values]) / len(xyz_values)

  xx = sum([(d[0] - x0) ** 2 for d in xyz_values])
  xy = sum([(d[0] - x0) * (d[1] - y0) for d in xyz_values])
  yy = sum([(d[1] - y0) ** 2 for d in xyz_values])
  xz = sum([(d[0] - x0) * (d[2] - z0) for d in xyz_values])
  yz = sum([(d[1] - y0) * (d[2] - z0) for d in xyz_values])

  A = matrix.sqr((xx, xy, 0.0, xy, yy, 0.0, 0.0, 0.0, len(xyz_values)))
  b = matrix.col((xz, yz, 0.0))
  x = A.inverse() * b
  c = z0 - x[0] * x0 - x[1] * y0
  return x[0], x[1], c

def smooth_invert_calibration_table(table, direction):
  '''Make smoothed version of the incoming table using plane linear interpolation
  and return the smoothed table and it's inverse. Scale for smoothing is relatively
  arbitrary.'''

  assert(direction in ['fast', 'slow'])

  smooth_scale = 3

  import copy

  # first duplicate the table

  smooth_table = copy.deepcopy(table)
  smooth_table_inverse = copy.deepcopy(table)

  # now smooth it out nicely

  size = table.focus()

  for j in range(size[0]):
    print j
    for i in range(size[1]):
      xyz = []
      for _j in range(j - smooth_scale, j + smooth_scale + 1):
        if _j < 0:
          continue
        if _j >= size[0]:
          continue
        for _i in range(i - smooth_scale, i + smooth_scale + 1):
          if _i < 0:
            continue
          if _i >= size[1]:
            continue
          xyz.append((_i, _j, table[_j, _i]))
      a, b, c = fit_plane(xyz)
      o = a * i + b * j + c
      smooth_table[j, i] = o
      if direction == 'fast':
        smooth_table_inverse[j, i] = - (a * (i + o) + b * j + c)
      else:
        smooth_table_inverse[j, i] = - (a * i + b * (j + o) + c)

  return smooth_table, smooth_table_inverse

xds_template = '''JOB=XYCORR
DETECTOR=PILATUS MINIMUM_VALID_PIXEL_VALUE=0 OVERLOAD=%(overload)d
DIRECTION_OF_DETECTOR_X-AXIS=%(fast_x).3f %(fast_y).3f %(fast_z).3f
DIRECTION_OF_DETECTOR_Y-AXIS=%(slow_x).3f %(slow_y).3f %(slow_z).3f
TRUSTED_REGION=0.0 1.41
NX=%(n_fast)d NY=%(n_slow)d QX=%(pixel_fast).4f QY=%(pixel_slow).4f
DETECTOR_DISTANCE=%(distance).2f
X-RAY_WAVELENGTH=%(wavelength).6f
INCIDENT_BEAM_DIRECTION=%(beam_x).3f %(beam_y).3f %(beam_z).3f
SENSOR_THICKNESS= %(thickness).3f
ORGX=%(origin_fast).2f ORGY=%(origin_slow).2f'''

def image_to_XDS_XYCORR(image_filename, sensor_thickness_mm):
  '''Generate an XYCORR input file from an image header via dxtbx, noting
  well that this will *tell lies* as the image is rescaled to give a 1:1
  correction table in 0.025 rather than 0.1 (original) pixel increments.'''

  from dxtbx import load
  from scitbx import matrix

  image = load(image_filename)

  beam = matrix.col(image.get_beam().get_s0())
  wavelength = image.get_beam().get_wavelength()
  d = image.get_detector()[0]
  fast = matrix.col(d.get_fast_axis())
  slow = matrix.col(d.get_slow_axis())
  normal = matrix.col(d.get_normal())
  origin = matrix.col(d.get_origin())
  distance = origin.dot(normal)
  offset = distance * normal - origin
  offset_fast = offset.dot(fast)
  offset_slow = offset.dot(slow)

  trusted = d.get_trusted_range()

  pixel_size = d.get_pixel_size()

  # this is in order slow, fast i.e. C order
  image_size = image.get_raw_data().focus()

  open('XDS.INP', 'w').write(xds_template % {
    'overload':trusted[1],
    'fast_x':fast.elems[0],
    'fast_y':fast.elems[1],
    'fast_z':fast.elems[2],
    'slow_x':slow.elems[0],
    'slow_y':slow.elems[1],
    'slow_z':slow.elems[2],
    'n_fast':image_size[1] * 4,
    'n_slow':image_size[0] * 4,
    'pixel_fast':pixel_size[0] / 4.0,
    'pixel_slow':pixel_size[1] / 4.0,
    'distance':distance,
    'wavelength':wavelength,
    'beam_x':beam.elems[0],
    'beam_y':beam.elems[1],
    'beam_z':beam.elems[2],
    'thickness':sensor_thickness_mm,
    'origin_fast':offset_fast / (pixel_size[0] / 4.0),
    'origin_slow':offset_slow / (pixel_size[1] / 4.0)
    })

  output = run_job('xds_par')

  # now read the correction tables in and scale back to pixels

  x_corrections_parallax = read_xds_calibration_file(
    'X-CORRECTIONS.cbf').as_double() / 40.0
  y_corrections_parallax = read_xds_calibration_file(
    'Y-CORRECTIONS.cbf').as_double() / 40.0

  # smooth and invert them...

  import matplotlib
  matplotlib.use('Agg')
  from matplotlib import pyplot
  pyplot.imshow(x_corrections_parallax.as_numpy_array())
  pyplot.savefig('x-corrections.png')
  pyplot.imshow(y_corrections_parallax.as_numpy_array())
  pyplot.savefig('y-corrections.png')

  print 'original maps made'

  x_map, x_inv = smooth_invert_calibration_table(x_corrections_parallax, 'fast')
  pyplot.imshow(x_map.as_numpy_array())
  pyplot.savefig('x-corrections-smooth.png')
  pyplot.imshow(x_inv.as_numpy_array())
  pyplot.savefig('x-corrections-invert.png')
  y_map, y_inv = smooth_invert_calibration_table(y_corrections_parallax, 'slow')
  pyplot.imshow(y_map.as_numpy_array())
  pyplot.savefig('y-corrections-smooth.png')
  pyplot.imshow(y_inv.as_numpy_array())
  pyplot.savefig('y-corrections-invert.png')



  # now read in the DECTRIS correction tables for this detector (assuming that it
  # is serial No. 100 P6M at Diamond...)




if __name__ == '__main__':
  import sys
  image_to_XDS_XYCORR(sys.argv[1], float(sys.argv[2]))
