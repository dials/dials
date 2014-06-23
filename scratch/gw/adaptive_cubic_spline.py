class bicubic_bilinear_spline(object):
  '''A class to provide support for bicubic spline interpolation falling back
  on bilinear spline interpolation at the corners. To avoid the fallback, pad
  the incoming parameter array appropriately.'''

  def __init__(self, point_array):
    self._point_array = point_array
    focus = point_array.focus()
    self._ny, self._nx = focus[0] - 1, focus[1] - 1

    from scitbx.array_family import flex
    self._44 = flex.grid(4, 4)
    self._aij_pij = self._make_aij_pij_coefficients_table()
    return

  def _make_aij_pij_coefficients_table(self):
    from scitbx.array_family import flex
    aij_pij = flex.double((
      (0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00,
       0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
      (0.00, 0.00, 0.00, 0.00, -0.50, 0.00, 0.50, 0.00,
       0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
      (0.00, 0.00, 0.00, 0.00, 1.00, -2.50, 2.00, -0.50,
       0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
      (0.00, 0.00, 0.00, 0.00, -0.50, 1.50, -1.50, 0.50,
       0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
      (0.00, -0.50, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
       0.00, 0.50, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
      (0.25, 0.00, -0.25, 0.00, 0.00, 0.00, 0.00, 0.00,
       -0.25, 0.00, 0.25, 0.00, 0.00, 0.00, 0.00, 0.00),
      (-0.50, 1.25, -1.00, 0.25, 0.00, 0.00, 0.00, 0.00,
       0.50, -1.25, 1.00, -0.25, 0.00, 0.00, 0.00, 0.00),
      (0.25, -0.75, 0.75, -0.25, 0.00, 0.00, 0.00, 0.00,
       -0.25, 0.75, -0.75, 0.25, 0.00, 0.00, 0.00, 0.00),
      (0.00, 1.00, 0.00, 0.00, 0.00, -2.50, 0.00, 0.00,
       0.00, 2.00, 0.00, 0.00, 0.00, -0.50, 0.00, 0.00),
      (-0.50, 0.00, 0.50, 0.00, 1.25, 0.00, -1.25, 0.00,
       -1.00, 0.00, 1.00, 0.00, 0.25, 0.00, -0.25, 0.00),
      (1.00, -2.50, 2.00, -0.50, -2.50, 6.25, -5.00, 1.25,
       2.00, -5.00, 4.00, -1.00, -0.50, 1.25, -1.00, 0.25),
      (-0.50, 1.50, -1.50, 0.50, 1.25, -3.75, 3.75, -1.25,
       -1.00, 3.00, -3.00, 1.00, 0.25, -0.75, 0.75, -0.25),
      (0.00, -0.50, 0.00, 0.00, 0.00, 1.50, 0.00, 0.00,
       0.00, -1.50, 0.00, 0.00, 0.00, 0.50, 0.00, 0.00),
      (0.25, 0.00, -0.25, 0.00, -0.75, 0.00, 0.75, 0.00,
       0.75, 0.00, -0.75, 0.00, -0.25, 0.00, 0.25, 0.00),
      (-0.50, 1.25, -1.00, 0.25, 1.50, -3.75, 3.00, -0.75,
       -1.50, 3.75, -3.00, 0.75, 0.50, -1.25, 1.00, -0.25),
      (0.25, -0.75, 0.75, -0.25, -0.75, 2.25, -2.25, 0.75,
       0.75, -2.25, 2.25, -0.75, -0.25, 0.75, -0.75, 0.25)))

    return aij_pij

  def _aij(self, pij):
    '''Derive coefficient matrix aij from 4 x 4 matrix pij'''

    assert(pij.focus() == (4, 4))

    # compute aij from pij by matrix operations on 16-vectors not 4x4 matrices

    aij = self._aij_pij.matrix_multiply((pij.as_1d()))

    return aij

  def _evaluate_bilinear(self, x, y):
    x0 = int(x)
    y0 = int(y)
    _x = x - x0
    _y = y - y0

    f00 = self._point_array[y0, x0]
    f01 = self._point_array[y0 + 1, x0]
    f10 = self._point_array[y0, x0 + 1]
    f11 = self._point_array[y0 + 1, x0 + 1]

    return f00 + (f10 - f00) * _x + (f01 - f00) * _y + \
      (f00 - f10 - f01 + f11) * _x * _y

  def _evaluate_bicubic(self, x, y):
    x0 = int(x)
    y0 = int(y)
    _x = x - x0
    _y = y - y0

    from scitbx.array_family import flex
    X = flex.double(flex.grid(4, 4))

    for i in range(4):
      for j in range(4):
        X[i,j] = (_x ** j) * (_y ** j)

    aij = self._aij(self._point_array[y0-1:y0+3,x0-1:x0+3])

    return aij.dot(X.as_1d())

  def evaluate(self, x, y):
    assert(x > 0 and x < self._nx)
    assert(y > 0 and y < self._ny)

    if x > 1 and x < (self._nx - 1) and \
      y > 1 and y < (self._ny - 1):
      return self._evaluate_bicubic(x, y)

    return self._evaluate_bilinear(x, y)

  def __call__(self, x, y):
    return self.evaluate(x, y)

def tst_bicubic_binear_spline(n_points=100):
  from scitbx.array_family import flex
  import math

  scale = 1.0 / n_points

  cos_x_sin_y = flex.double(flex.grid(n_points + 1, n_points + 1))
  for i in range(n_points + 1):
    for j in range(n_points + 1):
      cos_x_sin_y[i,j] = math.sin(math.pi * i * scale) * \
        math.cos(math.pi * j * scale)

  bbs = bicubic_bilinear_spline(cos_x_sin_y)

  import random

  n = 1000
  s = 0.0
  for j in range(n):
    x = n_points * random.random()
    y = n_points * random.random()
    cxsy = math.sin(math.pi * y * scale) * math.cos(math.pi * x * scale)
    s += (bbs(x, y) - cxsy) ** 2

  return math.sqrt(s / n)

if __name__ == '__main__':
  for n in 10, 20, 40, 100, 200, 400:
    print n, tst_bicubic_binear_spline(n)
