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

  derive a smoothed atenuation coefficient at a given energy in KeV, in cm ^ -1'''

  coefficients = [(3.0, 2217.228), (4.0, 1031.491), (5.0, 559.2),
                  (6.0, 335.287), (8.0, 147.0929), (10.0, 76.6337),
                  (15.0, 22.82002), (20.0, 9.49708)]

  assert(energy_kev >= 3.0)
  assert(energy_kev <= 20.0)

  for j, e_mu in enumerate(coefficients):
    e, mu = e_mu
    if e >= energy_kev:
      e_mu0 = coefficients[j-1]
      e_mu1 = coefficients[j]
      return log_interpolate(e_mu0[0], e_mu0[1], e_mu1[0], e_mu1[1], energy_kev)

  raise RuntimeError, 'cannot reach this point'

def compute_offset(t0, theta, mu):
  import math
  t = t0 / math.cos(theta)
  offset = math.sin(theta) * (1 - (1 + mu * t) * math.exp(- mu * t)) / \
    (mu * (1 - math.exp(- mu * t)))
  return offset

def dqe(t0, theta, mu):
  '''Compute DQE for the given thickness of sensor, for the given angle and linear
  absorption coefficient.'''

  import math

  t = t0 / math.cos(theta)

  return 1.0 - math.exp(-mu * t)

def work():
  '''320 micron sensor, 12.7 KeV photons, theta intersection angle 0 to 45 degrees
  by way of a code test, return values in pixels i.e. multiples of 172 microns.'''

  import math

  # all calculations performed in cm

  d2r = math.pi / 180.0
  mu_cm_127 = derive_absorption_coefficient_Si(12.7)
  mu_cm_170 = derive_absorption_coefficient_Si(17.0)

  t0 = 0.032
  pixel = 0.0172

  for j in range(-45, 46):
    theta = d2r * j
    o_127 = compute_offset(t0, theta, mu_cm_127)
    o_170 = compute_offset(t0, theta, mu_cm_170)
    print j, o_127 / pixel, o_170 / pixel

def work_dqe():
  import math
  t0 = 0.032
  for j in range(30, 201, 1):
    energy_kev = 0.1 * j
    mu = derive_absorption_coefficient_Si(energy_kev)
    print energy_kev, dqe(t0, 0.0, mu), dqe(t0, math.pi / 4, mu)

def work_compare_2005_paper():
  '''Run with 12 KeV photons, 300 micron sensor, 217 micron pixel size to compare
  with 2005 Hulsen, Bronnimann & Eikenberry paper. OK with the current state of 
  the parallax corrections this now gives what look to be the same answers as 
  in this paper - pretty accurate but does not agree with XDS...'''

  import math

  # all calculations performed in cm

  d2r = math.pi / 180.0
  mu_cm = derive_absorption_coefficient_Si(12.0)

  t0 = 0.03
  pixel = 0.0217

  for j in range(0, 61):
    theta = d2r * j
    o = compute_offset(t0, theta, mu_cm)
    print j, o / pixel

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

  assert(length == fast * slow)

  pixel_values = uncompress(packed = data[data_offset:],
                            fast = fast, slow = slow)

  return pixel_values

def validate_against_xds(xds_directory):
  '''Will look in xds_directory for (i) XPARM.XDS and (ii) X- and Y-CORRECTIONS.cbf
  files, and will rederive from the former the values for the latter...'''

  from rstbx.cftbx.coordinate_frame_converter import coordinate_frame_converter
  import os
  import math

  # read and derive all of the model parameters

  xparm = os.path.join(xds_directory, 'XPARM.XDS')
  cfc = coordinate_frame_converter(xparm)
  fast = cfc.get_c('detector_fast')
  slow = cfc.get_c('detector_slow')
  origin = cfc.get_c('detector_origin')
  wavelength = cfc.get('wavelength')
  normal = fast.cross(slow)

  energy_kev = 12.3985 / wavelength
  mu = derive_absorption_coefficient_Si(energy_kev)

  detector_size = cfc.get('detector_size_fast_slow')
  pixel_size = cfc.get('detector_pixel_size_fast_slow')

  beam_centre_fast_slow = cfc.derive_beam_centre_pixels_fast_slow()
  beam_centre = beam_centre_fast_slow[0] * pixel_size[0] * fast + \
    beam_centre_fast_slow[1] * pixel_size[1] * slow + origin

  min_offset_x = 0
  max_offset_x = 0
  min_offset_y = 0
  max_offset_y = 0

  for k in range(0, detector_size[1], 20):
    for j in range(0, detector_size[0], 20):
      p = j * pixel_size[0] * fast + k * pixel_size[1] * slow + origin
      theta = p.angle(normal)

      # beware this piece of code is in cm not mm - blame particle physicists and
      # cgs units

      t0 = 0.032
      pixel = 0.0172

      offset = compute_offset(t0, theta, mu) / pixel

      # end weirdness

      offset_vector = (p - beam_centre).normalize() * offset
      offset_x = fast.dot(offset_vector)
      offset_y = slow.dot(offset_vector)
      if offset_x < min_offset_x:
        min_offset_x = offset_x
      if offset_x > max_offset_x:
        max_offset_x = offset_x
      if offset_y < min_offset_y:
        min_offset_y = offset_y
      if offset_y > max_offset_y:
        max_offset_y = offset_y

  print min_offset_x, max_offset_x
  print min_offset_y, max_offset_y

if __name__ == '__main__':
  import sys
  #validate_against_xds(sys.argv[1])
  #work_dqe()
  work_compare_2005_paper()
  #work()
