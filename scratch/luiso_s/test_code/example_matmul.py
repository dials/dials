from __future__ import division
from dials.scratch.luiso_s import tst_prod , write_2d
from scitbx.array_family import flex
a = flex.double(flex.grid(3, 3))
b = flex.double(flex.grid(3, 3))
import random
contr = 0
for xpos in range(3):
  for ypos in range(3):
    a[ypos, xpos] = random.random()
    b[ypos, xpos] = 1.0
    contr += 1.0

print "a ="
write_2d(a)

print "b ="
write_2d(b)




from scitbx import matrix
a_mat = a.as_scitbx_matrix()
b_mat = b.as_scitbx_matrix()

x_mat = a_mat * b_mat

x = x_mat.as_flex_double_matrix()

print "x = a * b ="
x_np = x.as_numpy_array()
print x_np
print "= as_flex"
write_2d(x)




y = tst_prod(a, b)

print "y = a * b ="
y_np = y.as_numpy_array()
print y_np
print "= as_flex"
write_2d(y)
