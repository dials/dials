from __future__ import division
to_do = '''

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> to do
A * x = B

[A] = matriz
[x] = vector columna
[B] = vextor columna
A ** (-1) * (A * x) =A ** (-1) * B

x = A ** (-1) * B


'''

from scitbx import matrix
from dials.scratch.luiso_s import tst_prod , write_2d
from scitbx.array_family import flex
import random
a = flex.double(flex.grid(3, 3))
b = flex.double(flex.grid(3, 1))


for ypos in range(3):
    b[ypos, 0] = random.random()
    for xpos in range(3):
        a[ypos, xpos] = random.random()


print "a ="
write_2d(a)

print "b ="
write_2d(b)

a_mat = a.as_scitbx_matrix()
b_mat = b.as_scitbx_matrix()

x_mat = a_mat.inverse() * b_mat

print "x = A ** (-1) * B"
x = x_mat.as_flex_double_matrix()
write_2d(x)
print "x(as flex)"
flex.show(x)
