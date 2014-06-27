from __future__ import division

from mu_Si import derive_absorption_coefficient_Si
from parallax_xds import generate_xds_corrections

def compute_absolute_offset(t0, theta, mu):
  '''N.B. to use needs geometric term e.g. sin(theta)'''
  import math
  t = t0 / math.cos(theta)
  _mu = 1.0 / mu
  return (_mu - (t + _mu) * math.exp(- mu * t))

def compute_absolute_offset_xds(t0, theta, mu):
  import math
  return min(t0 / math.cos(theta), 1.0 / mu)

def generate_dials_corrections(image_filename, sensor_thickness_mm, 
                               energy_ev = None, method=compute_absolute_offset):
  '''Generate equivalent correction tables equivalent to those from XDS, but using
  the equations above.'''

  from dxtbx import load
  from scitbx import matrix
  from scitbx.array_family import flex
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

  fast_parallax = flex.double(flex.grid(image_size))
  slow_parallax = flex.double(flex.grid(image_size))
  
  for i in range(image_size[0]):
    for j in range(image_size[1]):
      p = (origin + i * S + j * F).normalize()
      theta = p.angle(normal)
      dot_f = - p.dot(fast)
      dot_s = - p.dot(slow)
      offset = method(sensor_thickness_mm, theta, mu)
      fast_parallax[i, j] = dot_f * offset
      slow_parallax[i, j] = dot_s * offset

  return fast_parallax, slow_parallax

def work(image_filename, sensor_thickness_mm, energy_ev=None):
  '''Exercise the DIALS implementation of the XDS correction, compare with 
  the XDS correction tables.'''
  from parallax_xds import generate_xds_corrections
  xds_parallax_x, xds_parallax_y = generate_xds_corrections(
    image_filename, sensor_thickness_mm, energy_ev)
  print min(xds_parallax_x), min(xds_parallax_y), \
    max(xds_parallax_x), max(xds_parallax_y)
  dials_parallax_x, dials_parallax_y = generate_dials_corrections(
    image_filename, sensor_thickness_mm, energy_ev, 
    method=compute_absolute_offset_xds)
  print min(dials_parallax_x), min(dials_parallax_y), \
    max(dials_parallax_x), max(dials_parallax_y)
  dx = xds_parallax_x - dials_parallax_x
  dy = xds_parallax_y - dials_parallax_y

  print min(dx), min(dy), max(dx), max(dy)
  
  import matplotlib
  matplotlib.use('Agg')
  from matplotlib import pyplot
  pyplot.imshow(dx.as_numpy_array())
  pyplot.colorbar()
  pyplot.savefig('dx.png')
  pyplot.imshow(dy.as_numpy_array())
  pyplot.savefig('dy.png')
  
  return

if __name__ == '__main__':
  import sys
  if len(sys.argv) == 3:
    work(sys.argv[1], float(sys.argv[2]))
  else:
    work(sys.argv[1], float(sys.argv[2]), energy_ev = int(sys.argv[3]))
    
