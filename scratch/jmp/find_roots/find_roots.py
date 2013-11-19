from __future__ import division


def bisection(func, x0, x1, eps, max_iter=50):

  # Evaluate function at brackets
  f_x0 = func(x0)
  f_x1 = func(x1)

  # Ensure function crosses zero
  if f_x0 * f_x1 >= 0.0:
    raise RuntimeError('Root is not bracketed')

  # Setup the difference and starting value
  if f_x0 < 0.0:
    dx = x1 - x0
    x = x0
  else:
    dx = x0 - x1
    x = x1

  # Iterate the given number of times
  for i in range(max_iter):
    dx *= 0.5
    xm = x + dx
    f_xm = func(xm)
    if f_xm <= 0.0:
      x = xm
    if abs(dx) < eps or f_xm == 0.0:
      return x

  # Raise an exception if we reach maximum iterations
  raise RuntimeError('Maximum iterations reached')


def secant(func, x0, x1, eps, max_iter=50):

  # Evaluate the function at the interval
  f_l = func(x0)
  f_x = func(x1)

  # Set the initial x
  if abs(f_l) < abs(f_x):
    x = x0
    xl = x1
    f_l, f_x = f_x, f_l
  else:
    xl = x0
    x = x1

  # Iterate
  for i in range(max_iter):
    dx = (xl - x) * f_x / (f_x - f_l)
    xl = x
    f_l = f_x
    x += dx
    f_x = func(x)

    # If we've reached the root, return
    if abs(dx) < eps or f_x == 0.0:
      return x

  # Raise an exception if we reach maximum iterations
  raise RuntimeError('Maximum iterations reached')


def regula_falsi(func, x0, x1, eps, max_iter=50):

  fl = func(x0)
  fh = func(x1)

  if fl * fh > 0.0:
    raise RuntimeError('Root is not bracketed')

  if fl < 0.0:
    xl = x0
    xh = x1
  else:
    xl = x1
    xh = x0
    fl, fh = fh, fl

  dx = xh - xl

  for i in range(max_iter):
    rtf = xl + dx * fl / (fl - fh)
    f = func(rtf)
    if f < 0.0:
      de = xl - rtf
      xl = rtf
      fl = f
    else:
      de = xh - rtf
      xh = rtf
      fh = f

    dx = xh - xl
    if abs(de) < eps or f == 0.0:
      return rtf

  # Raise an exception if we reach maximum iterations
  raise RuntimeError('Maximum iterations reached')


def sign(a, b):
  s = 1
  if b < 0.0:
    s = -1
  return s * abs(a)


def ridder(func, x0, x1, eps, max_iter=50):
  from math import sqrt

  # Evaluate function at interval
  fl = func(x0)
  fh = func(x1)

  # Ensure root is brackted
  if fl * fh > 0.0:
    raise RuntimeError('Root is not bracketed')
  elif fl == 0.0:
    return x0
  elif fh == 0.0:
    return x1

  xl = x0
  xh = x1
  ans = -9.99e99

  for i in range(max_iter):
    xm = 0.5 * (xl + xh)
    fm = func(xm)
    s = sqrt(fm * fm - fl * fh)
    if fl >= fh:
      sn = 1.0
    else:
      sn = -1.0
    xnew = xm + (xm - xl) * sn * fm / s
    if abs(xnew - ans) <= eps:
      return ans
    ans = xnew
    fnew = func(ans)
    if fnew == 0.0:
      return ans
    if sign(fm, fnew) != fm:
      xl = xm
      fl = fm
      xh = ans
      fh = fnew
    elif sign(fl, fnew) != fl:
      xh = ans
      fh = fnew
    elif sign(fh, fnew) != fh:
      xl = ans
      fl = fnew
    else:
      raise RuntimeError('How\'d I get here!')

    if abs(xh - xl) < eps:
      return ans

  # Raise an exception if we reach maximum iterations
  raise RuntimeError('Maximum iterations reached')


def brent(func, x1, x2, eps, max_iter=50):

  # Set the initial values
  a = x1
  b = x2
  c = x2
  fa = func(a)
  fb = func(b)
  fc = fb

  # Ensure root is brackted
  if fa * fb > 0.0:
    raise RuntimeError('Root is not bracketed')


  for i in range(max_iter):

    if (fb > 0.0 and fc > 0.0) or (fb < 0.0 and fc < 0.0):
      c = a
      fc = fa
      d = b - a
      e = b - a

    if abs(fc) < abs(fb):
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa

    tol1 = 2.0 * eps * abs(b) + 0.5 * eps
    xm = 0.5 * (c - b)

    if abs(xm) <= tol1 or fb == 0.0:
      return b
    if abs(e) >= tol1 and abs(fa) > abs(fb):
      s = fb / fa
      if a == c:
        p = 2.0 * xm * s
        q = 1.0 - s
      else:
        q = fa / fc
        r = fb / fc
        p = s * (2.0 * xm * q * (q - r) - (b - a)*(r - 1.0))
        q = (q - 1.0) * (r - 1.0)*(s - 1.0)

      if p > 0.0:
        q = -q

      p = abs(p)

      min1 = 3.0 * xm * q - abs(tol1 * q)
      min2 = abs(e * q)
      min3 = min([min1, min2])

      if 2.0 * p < min3:
        e = d
        d = p / q
      else:
        d = xm
        e = d

    else:
      d = xm
      e = d

    a = b
    fa = fb

    if (abs(d) > tol1):
      b += d
      fb = func(b)
    else:
      b += sign(tol1, xm)
      fb = func(b)

  # Raise an exception if we reach maximum iterations
  raise RuntimeError('Maximum iterations reached')


def newton(func, dfunc, x, eps, max_iter=50):

  # Iterate
  for i in range(max_iter):

    # Evaluate the function and its derivative at the
    # current value of x
    f_x = func(x)
    df_x = dfunc(x)

    # Ensure the derivative is non-zero
    if df_x == 0.0:
      raise ValueError('Derivative equals zero')

    # Calculate the change in x and x
    dx = f_x / df_x
    x -= dx

    # If convergence has been achieved then return
    if abs(dx) < eps:
      return x

  # Raise an exception if we reach maximum iterations
  raise RuntimeError('Maximum iterations reached')


def func(x):
  from math import sqrt
  d = 10
  l = 200
  h = 200
  return d - x * (1 + l / sqrt(h**2 + x**2))

def dfunc(x):
  from math import sqrt
  d = 10
  l = 200
  h = 200
  return -h*h*l / sqrt(h**2 + x**2)**3 - 1

def time_func(func):
  from time import time
  st = time()
  for i in range(10000):
    v = func()
  print v, time() - st

time_func(lambda: bisection(func, 0, 200, 1e-10))
time_func(lambda: secant(func, 0, 200, 1e-10))
time_func(lambda: regula_falsi(func, 0, 200, 1e-10))
time_func(lambda: ridder(func, 0, 200, 1e-10))
time_func(lambda: brent(func, 0, 200, 1e-10))
time_func(lambda: newton(func, dfunc, 0, 1e-10))
