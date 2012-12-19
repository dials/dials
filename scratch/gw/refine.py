from __future__ import division
import numpy
from scipy.optimize import minimize

def f(xy):
    from math import sin, cos, pow, pi
    x, y = xy
    return pow(x * x + y * y + pi * sin(x) + pi * cos(y), 2)

def J(xy):
    from math import sin, cos, pow, pi
    x, y = xy
    return numpy.array(
        [(4*x + 2*pi*cos(x))*(x**2 + y**2 + pi*sin(x) + pi*cos(y)),
         (4*y - 2*pi*sin(y))*(x**2 + y**2 + pi*sin(x) + pi*cos(y))])

def H(xy):
    from math import sin, cos, pow, pi
    x, y = xy
    return numpy.matrix(
        [[(2*x + pi*cos(x))*(4*x + 2*pi*cos(x)) +
          (-2*pi*sin(x) + 4)*(x**2 + y**2 + pi*sin(x) + pi*cos(y)),
            (4*x + 2*pi*cos(x))*(2*y - pi*sin(y))],
            [(2*x + pi*cos(x))*(4*y - 2*pi*sin(y)),
             (2*y - pi*sin(y))*(4*y - 2*pi*sin(y)) +
                (-2*pi*cos(y) + 4)*(x**2 + y**2 + pi*sin(x) + pi*cos(y))]])

def run_nelder_mead():
    xy0 = numpy.array([1.0, 1.0])
    res = minimize(f, xy0, method = 'nelder-mead')
    x, y = res.x
    return x, y, f((x, y))

def run_bfgs():
    xy0 = numpy.array([1.0, 1.0])
    res = minimize(f, xy0, method = 'BFGS', jac = J)
    x, y = res.x
    return x, y, f((x, y))

def run_newton_cg():
    xy0 = numpy.array([1.0, 1.0])
    res = minimize(f, xy0, method = 'Newton-CG', jac = J, hess = H)
    x, y = res.x
    return x, y, f((x, y))

import time
t0 = time.time()
for j in range(1000): run_nelder_mead()
t1 = time.time()
for j in range(1000): run_bfgs()
t2 = time.time()
for j in range(1000): run_newton_cg()
t3 = time.time()

print 'Nelder-Mead; %.2f' % (t1 - t0)
print 'BFGS:        %.2f' % (t2 - t1)
print 'Newton CG:   %.2f' % (t3 - t2)
