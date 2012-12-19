from __future__ import division
import sympy

symbols = 'x y z'
functions = 'sin cos pi exp'
f = '(z ** 2 + x ** 2 + y ** 2 + sin(x) + cos(y)) ** 2'

_syms = sympy.symbols(symbols)

fout = open('fJH.py', 'w')
fout.write('def f(%s):\n' % ''.join(symbols.split()))
fout.write('    import numpy\n')
fout.write('    from math import %s\n' % ', '.join(functions.split()))
fout.write('    %s = %s\n' % (', '.join(symbols.split()),
                              ''.join(symbols.split())))
fout.write('    return %s\n' % f)
fout.write('def J(%s):\n' % ''.join(symbols.split()))
fout.write('    import numpy\n')
fout.write('    from math import %s\n' % ', '.join(functions.split()))
fout.write('    %s = %s\n' % (', '.join(symbols.split()),
                              ''.join(symbols.split())))
fout.write('    return numpy.array([%s])\n' % \
           ', '.join([str(sympy.diff(f, _s)) for _s in _syms]))
fout.write('def H(%s):\n' % ''.join(symbols.split()))
fout.write('    import numpy\n')
fout.write('    from math import %s\n' % ', '.join(functions.split()))
fout.write('    %s = %s\n' % (', '.join(symbols.split()),
                              ''.join(symbols.split())))
fout.write('    return numpy.matrix([%s])\n' % \
           (', '.join(
               ['[%s]' % ', '.join(
                   [str(sympy.diff(sympy.diff(f, _t), _s)) for _s in _syms])
                   for _t in _syms])))
fout.close()
