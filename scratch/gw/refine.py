from __future__ import division
import numpy
from scipy.optimize import minimize

from fJH import f, J, H

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

print run_nelder_mead()
print run_bfgs()
print run_newton_cg()
