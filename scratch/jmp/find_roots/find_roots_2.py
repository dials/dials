from __future__ import division

from scitbx import matrix

def forward(x, h, l):
  from math import sqrt
  return x + x * l / sqrt(h**2 + x**2)

def durand_kerner_iteration(poly, roots):
  for i in range(len(roots)):
    r_i = roots[i]
    f_r = poly(r_i)
    den = 1.0
    for j in range(len(roots)):
      if i != j:
        r_j = roots[j]
        den *= (r_i - r_j)
    r = r_i - f_r / den
    roots[i] = r
  return roots

def durand_kerner(poly, roots, eps=1e-10):

  import copy

  while True:
    old_roots = copy.deepcopy(roots)
    roots = durand_kerner_iteration(poly_func, roots)
    err = sum([abs(r1 - r2) for r1, r2 in zip(roots, old_roots)])
    if err <= eps:
      break

  return roots


def aberth_iteration(poly, deriv, roots):

  new_roots = roots
  for i in range(len(roots)):
    r_i = roots[i]
    f_r = poly(r_i)
    fp_r = deriv(r_i)

    num = f_r / fp_r

    s = 0.0
    for j in range(len(roots)):
      if i != j:
        r_j = roots[j]
        s += (1.0 / (r_i - r_j))
    den = 1.0 - num * s
    r = r_i -num / den
    new_roots[i] = r

  return new_roots


def aberth(poly, deriv, roots, eps=1e-10):

  import copy

  while True:
    old_roots = copy.deepcopy(roots)
    roots = aberth_iteration(poly_func, deriv, roots)
    err = sum([abs(r1 - r2) for r1, r2 in zip(roots, old_roots)])
    if err <= eps:
      break

  return roots


def laguerre_iteration(poly, deriv, deriv2, x, n, div_items):
  from cmath import sqrt

  div = 1.0
  for d in div_items:
    div *= (x - d)


  p_x = poly(x) / div
  pp_x = deriv(x) / div
  ppp_x =deriv2(x) / div
  G = pp_x / p_x
  H = G**2 - ppp_x / p_x

  den_n = sqrt((n - 1) * (n*H - G**2))
  if abs(G) < 0:
    den = G - den_n
  else:
    den = G + den_n

  a = n / den

  return x - a

def laguerre(poly, deriv, deriv2, x, n, eps=1e-10):

  import copy
  div = []
  x_first = x
  for i in range(n):
    x = x_first
    while True:
      x1 = laguerre_iteration(poly_func, deriv, deriv2, x, n-i, div)
      err = abs(x1 - x)
      if err <= eps:
        break
      x = x1

    print x1
    div.append(x1)
#        return x1


h = 200.0
d = 10.0
l = 0.1

p = -2.0*d
q = h**2 + d**2 -l**2
r = -2.0*d*h**2
s = (h*d)**2

print p, q, r, s

print forward(9.99500872462, h, l)


def poly_func(x):

  return x**4 + p*x**3 + q*x**2 + r*x + s

def poly_derivative(x):
  return 4*x**3 + 3*p*x**2 + 2*q*x + r

def poly_2nd_derivative(x):
  return 12*x**2 + 6*p*x + 2*q

roots = [(0.4 + 0.9j)**0,
         (0.4 + 0.9j)**1,
         (0.4 + 0.9j)**2,
         (0.4 + 0.9j)**3]

diff = []

from math import log
x = 10 + 10j

x = laguerre(poly_func, poly_derivative, poly_2nd_derivative, x, 4)

print x

#from matplotlib import pylab
#pylab.plot(diff)
#pylab.show()

#m = matrix.sqr((0, 0, 0, a0,
#                1, 0, 0, a1,
#                0, 1, 0, a2,
#                0, 0, 1, a3)).transpose()
#
#from scitbx.linalg import householder
#from scitbx.array_family import flex


#A = m

#for i in range(1000):

#    A = flex.double(A)
#    A.reshape(flex.grid(4, 4))
#
#
#    h = householder.householder_qr_decomposition(A, True)
#    Q = matrix.sqr(list(h.q()))
#    R = matrix.sqr(list(h.r))
#
#
#    Ai = Q.transpose() * matrix.sqr(list(A)) * Q
#
#    print Ai
#
#    A = Ai
