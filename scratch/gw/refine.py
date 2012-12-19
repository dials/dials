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
    res = minimize(f, xy0, method = 'nelder-mead', options = {
        'xtol':1.0e-8, 'disp':True})
    x, y = res.x
    print x, y, f((x, y))

def run_bfgs():
    xy0 = numpy.array([1.0, 1.0])
    res = minimize(f, xy0, method = 'BFGS', jac = J, options = {
        'disp':True})
    x, y = res.x
    print x, y, f((x, y))

def run_newton_cg():
    xy0 = numpy.array([1.0, 1.0])
    res = minimize(f, xy0, method = 'Newton-CG', jac = J, hess = H, options = {
        'disp':True})
    x, y = res.x
    print x, y, f((x, y))
    
run_nelder_mead()
run_bfgs()
run_newton_cg()
