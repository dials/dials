from __future__ import division
from dials.scratch.luiso_s import tst_prod , write_2d
from scitbx.array_family import flex
a = flex.double(flex.grid(3, 3))
b = flex.double(flex.grid(3, 3))
contr = 0
for npos in range(3):
    a[0, npos] = 1
    b[npos, 0] = 1

print "a ="
write_2d(a)

print "b ="
write_2d(b)

x = tst_prod(a, b)

print "x = a * b"
write_2d(x)

