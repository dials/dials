from __future__ import division
from dials.scratch.luiso_s import tst_ref_prod , write_2d
from scitbx.array_family import flex
a = flex.double(flex.grid(3, 3))
b = flex.double(flex.grid(3, 3))
import random
for xpos in range(3):
  for ypos in range(3):
    a[ypos, xpos] = random.random()
    b[ypos, xpos] = 1.0

print "a ="
write_2d(a)

print "b ="
write_2d(b)

c = tst_ref_prod(a, b)

print "c ="
write_2d(c)

print "______________________________________"


print "a ="
write_2d(a)

print "b ="
write_2d(b)

print "c ="
write_2d(c)
