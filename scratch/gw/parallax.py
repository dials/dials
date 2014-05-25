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

if __name__ == '__main__':
  work()
