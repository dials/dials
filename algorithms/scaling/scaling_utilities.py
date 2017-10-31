from dials.array_family import flex
import numpy as np

def calc_s2d(data_man, reflection_table):
  reflection_table['phi'] = (reflection_table['xyzobs.px.value'].parts()[2]
                             * data_man.experiments.scan.get_oscillation()[1])
  (s0x, s0y, s0z) = data_man.experiments.beam.get_s0()
  reflection_table['s2'] = reflection_table['s1'] - (s0x, s0y, s0z)
  reflection_table['s2d'] = apply_inverse_rotation(data_man, reflection_table['s2'],
                                                   reflection_table['phi'])

def apply_inverse_rotation(data_man, vec, phi):
  (r0, r1, r2) = vec.parts()
  (ux, uy, uz) = data_man.experiments.goniometer.get_rotation_axis()
  from math import pi
  c_ph = flex.double(np.cos(2.0*pi*phi/180.0))
  s_ph = flex.double(np.sin(2.0*pi*phi/180.0))
  rx = (((c_ph + ((ux**2) * (1.0 - c_ph))) * r0)
        + (((ux * uy * (1.0 - c_ph)) + (uz * s_ph)) * r1)
        + (((uz * ux * (1.0 - c_ph)) + (uy * s_ph)) * r2))
  ry = ((((ux * uy * (1.0 - c_ph)) - (uz * s_ph)) * r0)
        + ((c_ph + ((uy**2) * (1.0 - c_ph))) * r1)
        + (((uz * uy * (1.0 - c_ph)) + (ux * s_ph)) * r2))
  rz = ((((ux * uz * (1.0 - c_ph)) + (uy * s_ph)) * r0)
        + (((uy * uz * (1.0 - c_ph)) - (ux * s_ph)) * r1)
        + ((c_ph + ((uz**2) * (1.0 - c_ph))) * r2))
  s2d = zip(rx, ry, rz)
  return flex.vec3_double(s2d)

def sph_harm_table(reflection_table, lmax):
  from scitbx import math
  import math as pymath
  sph_harm_terms = flex.reflection_table()
  (x, y, z) = reflection_table['s2d'].parts()
  phi_list = flex.double(np.arctan2(y, x))
  theta_list = flex.double(np.arctan2((((x**2) + (y**2))**0.5), z))
  sqrt2 = pymath.sqrt(2)
  for l in range(1, lmax+1):
    lfg = math.log_factorial_generator(2 * l + 1)
    nsssphe = math.nss_spherical_harmonics(l, 50000, lfg)
    for m in range(-l, l+1):
      sph_harm_list = []
      for i, phi in enumerate(phi_list):
        theta = theta_list[i]
        Ylm = nsssphe.spherical_harmonic(l, abs(m), phi, theta)
        if m < 0:
          r = sqrt2 * ((-1) ** m) * Ylm.imag
        elif m == 0:
          assert Ylm.imag == 0.0
          r = Ylm.real
        else:
          r = sqrt2 * ((-1) ** m) * Ylm.real
        sph_harm_list.append(r)
      sph_harm_terms[str(l)+','+str(m)] = flex.double(sph_harm_list)
  return sph_harm_terms
