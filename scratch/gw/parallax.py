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
    o_127 = - (1.0 / mu_cm_127) * math.sin(theta) * \
      math.log(0.5 + 0.5 * math.exp(- mu_cm_127 * t0 / math.cos(theta)))
    o_170 = - (1.0 / mu_cm_170) * math.sin(theta) * \
      math.log(0.5 + 0.5 * math.exp(- mu_cm_170 * t0 / math.cos(theta)))
    print j, o_127 / pixel, o_170 / pixel

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

  t0 = 0.032
  pixel = 0.0172

  min_offset_x = 0
  max_offset_x = 0
  min_offset_y = 0
  max_offset_y = 0

  for k in range(0, detector_size[1], 10):
    for j in range(0, detector_size[0], 10):
      p = j * pixel_size[0] * fast + k * pixel_size[1] * slow + origin
      theta = p.angle(normal)
      offset = - (1.0 / mu) * math.sin(theta) * \
        math.log(0.5 + 0.5 * math.exp(- mu * t0 / math.cos(theta))) / pixel
      offset_vector = (p - beam_centre).normalize()
      offset_x = fast.dot(offset_vector * offset)
      offset_y = slow.dot(offset_vector * offset)
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
  validate_against_xds(sys.argv[1])
