from __future__ import division
def work_dqe():
  from parallax import dqe, derive_absorption_coefficient_Si
  t0 = 0.032
  for j in range(30, 201, 1):
    energy_kev = 0.1 * j
    mu = derive_absorption_coefficient_Si(energy_kev)
    print energy_kev, dqe(0.032, 0.0, mu), dqe(0.05, 0.0, mu), dqe(0.1, 0.0, mu)

def recover_xds_silicon(wavelength):
  '''SILICON=Fraction of intensity loss per mm due to absorption in silicon.
  In XDS - try to recover this value N.B. calculations still done in cm due to
  cgs nonsense...'''

  from parallax import derive_absorption_coefficient_Si
  import math
  energy_kev = 12.3985 / wavelength
  mu = derive_absorption_coefficient_Si(energy_kev)
  return 1 / (1 - mu * math.exp(- mu * 0.1))

if __name__ == '__main__':
  import sys
  print recover_xds_silicon(float(sys.argv[1]))
