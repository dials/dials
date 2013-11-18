from __future__ import division

def covar_contrib(x1, x2, weights):
  '''Compute contribution to covariance matrix.'''

  covar = sum([_x1 * _x2 * _w for _x1, _x2, _w in zip(x1, x2, weights)]) / \
      sum(weights)

  return covar

def form_covariance_matrix(pixel_list, _f, _r, _c):
  '''Compute covariance matrix.'''

  # form the variance matrix

  pixels = [(f - _f, r - _r, c - _c, d) for f, r, c, d in pixel_list]

  data = ([pixel[0] for pixel in pixels],
          [pixel[1] for pixel in pixels],
          [pixel[2] for pixel in pixels],
          [pixel[3] for pixel in pixels])

  m_elems = []

  for i in range(3):
    for j in range(3):
      m_elems.append(covar_contrib(data[i], data[j], data[3]))

  from scitbx.array_family import flex
  from scitbx.linalg import eigensystem
  from scitbx import matrix

  m = flex.double(flex.grid(3, 3))
  for j in range(9):
    m[j] = m_elems[j]

  s = eigensystem.real_symmetric(m)
  values = s.values()
  vector_elems = s.vectors()
  vectors = []

  for j in range(3):
    vectors.append(matrix.col(vector_elems[j * 3:j * 3 + 3]))

  return values, vectors

def nint(a):
  return int(round(a))

def generate_fake_profile(axes, variances, counts_total):
  '''Generate fake reflection profile with principle axes axes[0] ... [2],
  variance along each axis of variances[0..2] and a total of counts_total
  counts. Assume that the axes / variances are given in units of pixels.'''

  import math
  import random

  grid = { }

  sigmas = map(math.sqrt, variances)

  for j in range(counts_total):
    R = [random.gauss(0, s) for s in sigmas]
    x = R[0] * axes[0] + R[1] * axes[1] + R[2] * axes[2]
    c = tuple(map(nint, x))
    if not c in grid:
      grid[c] = 0
    grid[c] += 1

  pixels = []
  for c in grid:
    pixels.append((c[0], c[1], c[2], grid[c]))

  return pixels
