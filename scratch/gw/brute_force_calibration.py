from __future__ import division
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
DETECTOR=PILATUS MINIMUM_VALID_PIXEL_VALUE=0 OVERLOAD=1000000
DIRECTION_OF_DETECTOR_X-AXIS=1 0 0
DIRECTION_OF_DETECTOR_Y-AXIS=0 1 0
TRUSTED_REGION=0.0 1.41
NX=2463 NY=2527 QX=0.172 QY=0.172
DETECTOR_DISTANCE=%(distance).2f
X-RAY_WAVELENGTH=%(wavelength).6f
INCIDENT_BEAM_DIRECTION= 0 0 1
SENSOR_THICKNESS= %(thickness).3f
ORGX=1231.5 ORGY=1263.5
'''

def run_xds_xycorr(distance, wavelength, thickness):
  '''Run XDS XYCORR step for given distance, wavelength & thickness, read
  back the X, Y correction map, compute offset as a function of angle of
  incoming ray, return this.'''

  xds_output = xds_template % {'distance':distance,
                               'wavelength':wavelength,
                               'thickness':thickness}

  open('XDS.INP', 'w').write(xds_output)
  output = run_job('xds')
  return analyse_corrections(distance, wavelength, thickness)

def nint(a):
  return int(round(a))

def analyse_corrections(distance, wavelength, thickness):
  x_corr = read_xds_calibration_file('X-CORRECTIONS.cbf').as_double() * 0.1
  y_corr = read_xds_calibration_file('Y-CORRECTIONS.cbf').as_double() * 0.1
  o_squared = x_corr * x_corr + y_corr * y_corr

  dir_x = 2463 * 0.5
  dir_y = 2527 * 0.5

  pixel = 0.172

  from scitbx import matrix
  import math

  n = matrix.col((0, 0, 1))

  theta_o = { }

  for j in range(101):
    x = 0.01 * dir_x * j
    y = 0.01 * dir_y * j
    nx = nint(x / 4)
    ny = nint(y / 4)
    o = math.sqrt(o_squared[ny, nx])
    p = matrix.col((4 * pixel * (dir_x - x), 4 * pixel * (dir_y - y), distance))
    theta = p.angle(n)
    theta_o[theta] = o

  return theta_o

def main(distance=150):
  thickness = 0.32

  from parallax import derive_absorption_coefficient_Si, compute_offset, \
     compute_offset_dectris

  t0_cm = 0.032
  pixel_cm = 0.0172

  for energy_ev in range(8000, 16001, 4000):

    mu_cm = derive_absorption_coefficient_Si(energy_ev * 0.001)

    wavelength = 12.3985  / (energy_ev * 0.001)
    theta_o = run_xds_xycorr(distance, wavelength, thickness)

    fout = open('p%d.dat' % energy_ev, 'w')
    for theta in sorted(theta_o):
      model = compute_offset(t0_cm, theta, mu_cm) / pixel_cm
      model_d = compute_offset_dectris(t0_cm, theta, mu_cm) / pixel_cm
      fout.write('%.3f %.3f %.3f %.3f\n' % (theta, theta_o[theta], model, model_d))
    fout.close()

def main_plot(distance, energy_ev):
  thickness = 0.32

  from parallax import derive_absorption_coefficient_Si, compute_offset, \
     compute_offset_dectris

  t0_cm = 0.032
  pixel_cm = 0.0172

  mu_cm = derive_absorption_coefficient_Si(energy_ev * 0.001)

  wavelength = 12.3985  / (energy_ev * 0.001)
  theta_o = run_xds_xycorr(distance, wavelength, thickness)

  import math
  displacement = []
  xds = []
  m_gw = []
  m_dec = []

  for theta in sorted(theta_o):
    model = compute_offset(t0_cm, theta, mu_cm) / pixel_cm
    model_d = compute_offset_dectris(t0_cm, theta, mu_cm) / pixel_cm
    displacement.append(distance * math.tan(theta))
    xds.append(theta_o[theta])
    m_gw.append(model)
    m_dec.append(model_d)

  import matplotlib
  matplotlib.use('Agg')
  from matplotlib import pyplot
  p1, p2, p3 = pyplot.plot(displacement, xds, displacement, m_gw,
                           displacement, m_dec)
  pyplot.xlabel('Displacement, pixels')
  pyplot.ylabel('Offset, pixels')
  pyplot.title('Parallax correction offsets at %deV' % energy_ev)
  pyplot.legend([p1, p2, p3], ['Calculated by XDS', 'GW model', 'Dectris model'],
                loc=2)
  pyplot.savefig('%d.png' % energy_ev)

if __name__ == '__main__':
  import sys
  main_plot(150, energy_ev=float(sys.argv[1]))
