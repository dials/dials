import sympy
x, y = sympy.symbols('x y')
f = '(x*x + y*y + pi * sin(x) + pi * cos(y)) ** 2'
dfdx = sympy.diff(f, x)
dfdy = sympy.diff(f, y)
dfdxx = sympy.diff(dfdx, x)
dfdxy = sympy.diff(dfdx, y)
dfdyx = sympy.diff(dfdy, x)
dfdyy = sympy.diff(dfdy, y)

print 'Function'
print f
print 'Jacobian'
print 'return numpy.array([%s, %s])' % (dfdx, dfdy)
print 'Hessian'
print 'return numpy.matrix([[%s, %s], [%s, %s]])' % (dfdxx, dfdxy, dfdyx, dfdyy)
