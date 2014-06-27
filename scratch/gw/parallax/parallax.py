from __future__ import division
def log_interpolate(x0, y0, x1, y1, x):
  '''Return y(x) where we have fit a linear model to ln(y)(x)'''
  import math

  ly0 = math.log(y0)
  ly1 = math.log(y1)

  ly = ly0 + (ly1 - ly0) * (x - x0) / (x1 - x0)

  return math.exp(ly)

def derive_absorption_coefficient_Si(energy_kev):
  '''From data from

  http://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z14.html

  derive a smoothed atenuation coefficient at a given energy in KeV, in mm ^ -1'''

  if True:

    # computed from mu_en

    coefficients = [(3.0, 2217.228), (4.0, 1031.491), (5.0, 559.2),
                    (6.0, 335.287), (8.0, 147.0929), (10.0, 76.6337),
                    (15.0, 22.82002), (20.0, 9.49708)]

    assert(energy_kev >= 3.0)
    assert(energy_kev <= 20.0)

  else:

    # computed from mu

    coefficients = [(2.0, 6470.41), (3.0, 2279.672), (4.0, 1055.257),
                    (5.0, 570.85), (6.0, 342.51), (8.0, 150.7044),
                    (10.0, 78.9637), (15.0, 24.0922), (20.0, 10.40112),
                    (30.0, 3.34588), (40.0, 1.633796), (50.0, 1.021705),
                    (60.0, 0.747231)]

    assert(energy_kev >= 2.0)
    assert(energy_kev <= 60.0)

  for j, e_mu in enumerate(coefficients):
    e, mu = e_mu
    if e >= energy_kev:
      e_mu0 = coefficients[j-1]
      e_mu1 = coefficients[j]
      # / 10 to get per mm not per cm
      return 0.1 * log_interpolate(e_mu0[0], e_mu0[1], e_mu1[0], e_mu1[1],
                                   energy_kev)

  raise RuntimeError, 'cannot reach this point'

def compute_offset(t0, theta, mu):
  import math
  t = t0 / math.cos(theta)
  _mu = 1.0 / mu
  return math.sin(theta) * (_mu - (t + _mu) * math.exp(- mu * t))

def compute_offset_xds(t0, theta, mu):
  import math
  return min(t0 * math.tan(theta), math.sin(theta) / mu)

def dqe(t0, theta, mu):
  '''Compute DQE for the given thickness of sensor, for the given angle and linear
  absorption coefficient.'''

  import math
  t = t0 / math.cos(theta)
  return 1.0 - math.exp(-mu * t)

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

def image_to_XDS_XYCORR(image_filename, sensor_thickness_mm, energy_ev=None):
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

  # now read the correction tables in and scale back to pixels - recall pixel
  # size / 4 above..

  x_corrections_parallax = read_xds_calibration_file(
    'X-CORRECTIONS.cbf').as_double() / 40.0
  y_corrections_parallax = read_xds_calibration_file(
    'Y-CORRECTIONS.cbf').as_double() / 40.0

  if False:
    return x_corrections_parallax, y_corrections_parallax

  from scitbx.array_family import flex
  return flex.sqrt(x_corrections_parallax * x_corrections_parallax + \
                   y_corrections_parallax * y_corrections_parallax)

def image_to_test(image_filename, sensor_thickness_mm, energy_ev=None):
  xds_table = image_to_XDS_XYCORR(image_filename, sensor_thickness_mm,
                                  energy_ev=energy_ev)

  from dxtbx import load
  from scitbx import matrix
  import math
  import random
  import os

  image = load(image_filename)

  beam = matrix.col(image.get_beam().get_s0()).normalize()

  if energy_ev:
    energy_kev = energy_ev * 0.001
  else:
    wavelength = image.get_beam().get_wavelength()
    energy_kev = 12.3985 / wavelength

  mu = derive_absorption_coefficient_Si(energy_kev)
  d = image.get_detector()[0]
  fast = matrix.col(d.get_fast_axis())
  slow = matrix.col(d.get_slow_axis())
  normal = matrix.col(d.get_normal())
  origin = matrix.col(d.get_origin())
  distance = origin.dot(normal)
  offset = distance * normal - origin
  offset_fast = offset.dot(fast)
  offset_slow = offset.dot(slow)
  pixel_size = d.get_pixel_size()

  # this is in order slow, fast i.e. C order
  image_size = image.get_raw_data().focus()

  F = fast * pixel_size[0]
  S = slow * pixel_size[1]

  _theta = []
  _xds = []
  _dials = []
  _xds2 = []

  for j in range(10000):
    x = random.random() * image_size[1]
    y = random.random() * image_size[0]
    p = origin + y * S + x * F
    theta = p.angle(normal)
    xds = xds_table[int(y), int(x)]
    dials = compute_offset(sensor_thickness_mm, theta, mu)
    xds2 = compute_offset_xds(sensor_thickness_mm, theta, mu)

    _theta.append(theta * 180.0 / math.pi)
    _xds.append(xds)
    _dials.append(dials / pixel_size[0])
    _xds2.append(xds2 / pixel_size[0])

  # sort the results

  values = { }

  for t, x, d, x2 in zip(_theta, _xds, _dials, _xds2):
    values[t] = x, d, x2

  _theta, _xds, _dials, _xds2 = [], [], [], []
  for t in sorted(values):
    x, d, x2 = values[t]
    _theta.append(t)
    _xds.append(x)
    _dials.append(d)
    _xds2.append(x2)

  import matplotlib
  matplotlib.use('Agg')
  from matplotlib import pyplot
  p1, p2, p3 = pyplot.plot(_theta, _xds, _theta, _dials, _theta, _xds2)
  pyplot.xlabel('Angle, degrees')
  pyplot.ylabel('Offset, pixels')
  pyplot.title('Parallax correction offsets for %s' % os.path.split(
    image_filename)[-1])
  pyplot.legend([p1, p2, p3],
                ['Calculated by XDS', 'DIALS model', 'XDS model'],
                loc=2)
  if energy_ev:
    pyplot.savefig('corrections%d.png' % energy_ev)
  else:
    pyplot.savefig('corrections.png')

  return

if __name__ == '__make_xds_pictures__':
  import sys
  from scitbx.array_family import flex
  xds_parallax_x, xds_parallax_y = image_to_XDS_XYCORR(
    sys.argv[1], float(sys.argv[2]))

  import matplotlib
  matplotlib.use('Agg')
  from matplotlib import pyplot
  pyplot.imshow(xds_parallax_x.as_numpy_array())
  pyplot.colorbar()
  pyplot.savefig('xds_parallax_x.png')
  pyplot.imshow(xds_parallax_y.as_numpy_array())
  pyplot.savefig('xds_parallax_y.png')

if __name__ == '__main__':
  import sys
  if len(sys.argv) == 3:
    xds_parallax = image_to_test(sys.argv[1], float(sys.argv[2]))
  else:
    for arg in sys.argv[3:]:
      xds_parallax = image_to_test(sys.argv[1], float(sys.argv[2]),
                                   energy_ev = int(arg))
