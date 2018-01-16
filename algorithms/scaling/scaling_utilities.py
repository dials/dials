from dials.array_family import flex
import numpy as np
import cPickle as pickle

def calc_s2d(reflection_table, experiments):
  '''calculates diffraction vector in crystal frame'''
  reflection_table['phi'] = (reflection_table['xyzobs.px.value'].parts()[2]
                             * experiments.scan.get_oscillation()[1])
  (s0x, s0y, s0z) = experiments.beam.get_s0()
  reflection_table['s2'] = reflection_table['s1'] - (s0x, s0y, s0z)
  rot_axis = experiments.goniometer.get_rotation_axis()
  from math import pi
  angles = reflection_table['phi'] * -1.0 * pi / 180 #want to do an inverse rot.
  reflection_table['s2d'] = rotate_vectors_about_axis(rot_axis,
    reflection_table['s2'], angles)
  reflection_table['s2d'] = align_rotation_axis_along_z(rot_axis,
    reflection_table['s2d'])
  return reflection_table

def rotate_vectors_about_axis(rot_axis, vectors, angles):
  #assumes angles in radians
  (r0, r1, r2) = vectors.parts()
  (ux, uy, uz) = list(rot_axis)
  #normalise
  modulus = (ux**2 + uy**2 + uz**2)**0.5
  (ux, uy, uz) = (ux/modulus, uy/modulus, uz/modulus)
  from math import pi
  c_ph = flex.double(np.cos(angles))
  s_ph = flex.double(np.sin(angles))
  rx = (((c_ph + ((ux**2) * (1.0 - c_ph))) * r0)
        + (((ux * uy * (1.0 - c_ph)) - (uz * s_ph)) * r1)
        + (((uz * ux * (1.0 - c_ph)) + (uy * s_ph)) * r2))
  ry = ((((ux * uy * (1.0 - c_ph)) + (uz * s_ph)) * r0)
        + ((c_ph + ((uy**2) * (1.0 - c_ph))) * r1)
        + (((uz * uy * (1.0 - c_ph)) - (ux * s_ph)) * r2))
  rz = ((((ux * uz * (1.0 - c_ph)) - (uy * s_ph)) * r0)
        + (((uy * uz * (1.0 - c_ph)) + (ux * s_ph)) * r1)
        + ((c_ph + ((uz**2) * (1.0 - c_ph))) * r2))
  rotated_vectors = zip(rx, ry, rz)
  return flex.vec3_double(rotated_vectors)

def align_rotation_axis_along_z(exp_rot_axis, vectors):
  (ux, uy, uz) = list(exp_rot_axis)
  cross_prod_uz = (uy, -1.0*ux, 0.0)
  cpx, cpy, cpz = list(cross_prod_uz)
  from math import acos, pi
  angle_between_u_z = -1.0*acos(uz/((ux**2 + uy**2 + uz**2)**0.5))
  phi = flex.double([angle_between_u_z]*len(vectors))
  new_vectors = rotate_vectors_about_axis(cross_prod_uz, vectors, phi)
  return flex.vec3_double(new_vectors)

def sph_harm_table(reflection_table, lmax):
  from scitbx import math, sparse
  import math as pymath

  order = lmax
  lfg =  math.log_factorial_generator(2 * order + 1)
  n_params = 0
  for i in range(1,lmax+1):
    n_params += (2*i) +1
  #sph_harm_terms = flex.reflection_table()
  sph_harm_terms = sparse.matrix(len(reflection_table), n_params)
  (x, y, z) = reflection_table['s2d'].parts()

  phi_list = flex.double(np.arctan2(y, x))
  theta_list = flex.double(np.arctan2((((x**2) + (y**2))**0.5), z))
  #phi_list = flex.double(np.arctan2(z, y))
  #theta_list = flex.double(np.arctan2((((z**2) + (y**2))**0.5), x))
  sqrt2 = pymath.sqrt(2)
  nsssphe = math.nss_spherical_harmonics(order, 5000, lfg)
  counter = 0
  ziplist = zip(phi_list, theta_list)
  '''phi_round = phi_list.iround()
  theta_round = theta_list.iround()
  ziplist=zip(phi_round, theta_round)
  try:
    sphfile = open('sph_harm_table.pickle')
    shtable = pickle.load(sphfile)
    sphfile.close()
  except:
    generate_sph_harm_terms(4)
    sphfile = open('sph_harm_table.pickle')
    shtable = pickle.load(sphfile)
    sphfile.close()

  for i, (phi, theta) in enumerate(ziplist):
    j = theta + (phi*360)
    print(dir(shtable.rows()))
    print(list(shtable.rows()[0]))
    sph_harm_terms.assign_block(shtable[:][j], j, 0)'''


  for l in range(1, lmax+1):
    for m in range(-l, l+1):
      if m < 0:
        for i, (phi, theta) in enumerate(ziplist):
          sph_harm_terms[i, counter] = (sqrt2 * ((-1) ** m)
            * nsssphe.spherical_harmonic(l, -1*m, theta, phi).imag)
      elif m == 0:
        for i, (phi, theta) in enumerate(ziplist):
          sph_harm_terms[i, counter] = (
            nsssphe.spherical_harmonic(l, m, theta, phi).real)
      else:
        for i, (phi, theta) in enumerate(ziplist):
          sph_harm_terms[i, counter] = (sqrt2 * ((-1) ** m)
            * nsssphe.spherical_harmonic(l, m, theta, phi).real)
      counter += 1
  return sph_harm_terms


def generate_sph_harm_terms(lmax):
  from scitbx import math, sparse
  import math as pymath
  import numpy as np

  order = lmax
  lfg =  math.log_factorial_generator(2 * order + 1)
  n_params = 0
  for i in range(1,lmax+1):
    n_params += (2*i) +1
  sph_harm_terms = flex.reflection_table()
  for l in range(1, lmax+1):
    for m in range(-l, l+1):
      sph_harm_terms[str(l)+','+str(m)] = flex.double([0.0]*360*180)
  #sph_harm_test = sparse.matrix((360*180), n_params)
  sqrt2 = pymath.sqrt(2)
  nsssphe = math.nss_spherical_harmonics(order, 5000, lfg)
  counter = 0

  phis = flex.double(np.linspace(start=0, stop=180, num=180, endpoint=False)).iround()
  thetas = flex.double(np.linspace(start=0, stop=360, num=360, endpoint=False)).iround()
  print("generating spherical harmonic table")

  for l in range(1, lmax+1):
    for m in range(-l, l+1):
      col = str(l)+','+str(m)
      if m < 0:
        fac = float(((-1) ** m))
        for phi in phis:
          for theta in thetas:
            i = theta + (phi*360)
            sph_harm_terms[col][i] = sqrt2 * fac * nsssphe.spherical_harmonic(l, -1*m, theta, phi).imag
      elif m == 0:
        for phi in phis:
          for theta in thetas:
            i = theta + (phi*360)
            sph_harm_terms[col][i]  = nsssphe.spherical_harmonic(l, m, theta, phi).real
      else:
        fac = float(((-1) ** m))
        for phi in phis:
          for theta in thetas:
            i = theta + (phi*360)
            sph_harm_terms[col][i]  = (sqrt2 * fac
              * nsssphe.spherical_harmonic(l, m, theta, phi).real)
      counter += 1

  data_file = open('sph_harm_table.pickle', 'w')
  pickle.dump(sph_harm_terms, data_file)
  data_file.close()
