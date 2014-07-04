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

  # FIXME need to decide which absorption coefficient to use - I think the mu_en
  # one which refers to energy deposition. N.B. only tabulated likely range 3 -
  # 20 keV

  use_mu_en = True

  if use_mu_en:

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
      # * 0.1 to change from per cm to per mm
      return 0.1 * log_interpolate(e_mu0[0], e_mu0[1], e_mu1[0], e_mu1[1],
                                   energy_kev)

  raise RuntimeError, 'cannot reach this point'
