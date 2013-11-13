from __future__ import division
import numpy
from scipy.optimize import minimize

from fJH import f, J, H

def run_nelder_mead():
  xy0 = numpy.array([1.0, 1.0, 1.0])
  res = minimize(f, xy0, method = 'nelder-mead')
  x, y, z = res.x
  return x, y, z, f((x, y, z))

def run_bfgs():
  xy0 = numpy.array([1.0, 1.0, 1.0])
  res = minimize(f, xy0, method = 'BFGS', jac = J)
  x, y, z = res.x
  return x, y, z, f((x, y, z))

def run_newton_cg():
  xy0 = numpy.array([1.0, 1.0, 1.0])
  res = minimize(f, xy0, method = 'Newton-CG', jac = J, hess = H)
  x, y, z = res.x
  return x, y, z, f((x, y, z))

import time
t0 = time.time()
for j in range(1000): run_nelder_mead()
t1 = time.time()
for j in range(1000): run_bfgs()
t2 = time.time()
for j in range(1000): run_newton_cg()
t3 = time.time()

print 'Nelder-Mead: %2f' % (t1 - t0)
print 'BFGS:        %2f' % (t2 - t1)
print 'Newton-CG:   %2f' % (t3 - t2)
