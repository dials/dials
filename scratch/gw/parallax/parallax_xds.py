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

xds_template = '''JOB=XYCORR
DETECTOR=PILATUS MINIMUM_VALID_PIXEL_VALUE=0 OVERLOAD=%(overload)d
DIRECTION_OF_DETECTOR_X-AXIS=%(fast_x).3f %(fast_y).3f %(fast_z).3f
DIRECTION_OF_DETECTOR_Y-AXIS=%(slow_x).3f %(slow_y).3f %(slow_z).3f
TRUSTED_REGION=0.0 1.41
NX=%(n_fast)d NY=%(n_slow)d QX=%(pixel_fast).4f QY=%(pixel_slow).4f
DETECTOR_DISTANCE=%(distance).2f
X-RAY_WAVELENGTH=%(wavelength).6f
INCIDENT_BEAM_DIRECTION=%(beam_x).3f %(beam_y).3f %(beam_z).3f
SILICON= %(silicon).3f
SENSOR_THICKNESS= %(thickness).3f
ORGX=%(origin_fast).2f ORGY=%(origin_slow).2f
'''

# FIXME need export_xds_corrections()

def generate_xds_corrections(image_filename, sensor_thickness_mm,
                             energy_ev=None):
  '''Generate an XYCORR input file from an image header via dxtbx, noting
  well that this will *tell lies* as the image is rescaled to give a 1:1
  correction table in 0.025 rather than 0.1 (original) pixel increments.'''

  from dxtbx import load
  from scitbx import matrix

  image = load(image_filename)

  beam = matrix.col(image.get_beam().get_s0()).normalize()
  if energy_ev:
    wavelength = 12398.5 / energy_ev
  else:
    wavelength = image.get_beam().get_wavelength()
    energy_ev = 12398.5 / wavelength

  from mu_Si import derive_absorption_coefficient_Si
  silicon = derive_absorption_coefficient_Si(0.001 * energy_ev)
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
    'silicon':silicon,
    'thickness':sensor_thickness_mm,
    'origin_fast':offset_fast / (pixel_size[0] / 4.0),
    'origin_slow':offset_slow / (pixel_size[1] / 4.0)
    })

  output = run_job('xds_par')

  # now read the correction tables in and scale to mm - recall pixel size / 4
  # above..

  x_corrections_parallax = read_xds_calibration_file(
    'X-CORRECTIONS.cbf').as_double() * (pixel_size[1] / 40.0)
  y_corrections_parallax = read_xds_calibration_file(
    'Y-CORRECTIONS.cbf').as_double() * (pixel_size[0] / 40.0)

  return x_corrections_parallax, y_corrections_parallax
