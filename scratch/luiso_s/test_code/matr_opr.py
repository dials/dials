from __future__ import division
example = '''
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  from Graeme

>>> from scitbx import matrix
>>> import random
>>> m = matrix.sqr([random.random() for j in range(9)])
>>> m
matrix.rec(elems=(0.647908750543, ..., 0.874819390125), n=(3,3))
>>> m.inverse()
matrix.rec(elems=(4.57721132912, ..., 1.55733448049), n=(3,3))
>>> _ * m
matrix.rec(elems=(1.0, ..., 1.0), n=(3,3))
>>> i = m.inverse()
>>> m * i
matrix.rec(elems=(1.0, ..., 1.0), n=(3,3))
>>> i * m
matrix.rec(elems=(1.0, ..., 1.0), n=(3,3))
>>>


>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   from me



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


>>>>>>>>>>>>>>>>>>>>>>>>. from James


from scitbx import matrix

a = matrix.sqr((0, 0, 0, 0))

a.as_flex_double_matrix()
Out[47]: <scitbx_array_family_flex_ext.double at 0x455bfc8>

b = a.as_flex_double_matrix()

b.as_scitbx_matrix()
Out[50]: matrix.rec(elems=(0.0, ..., 0.0), n=(2,2))

'''

from scitbx import matrix
from dials.scratch.luiso_s import tst_prod , write_2d
from scitbx.array_family import flex
a = flex.double(flex.grid(3, 3))
#b = flex.double(flex.grid(3, 3))
import random
for xpos in range(3):
    for ypos in range(3):
        a[ypos, xpos] = random.random()
#        b[ypos, xpos] = 1.0

print "a ="
write_2d(a)

#print "b ="
#write_2d(b)

a_mat = a.as_scitbx_matrix()

print a_mat.inverse()
x_mat = a_mat.inverse()
print "inverse(a)= "
x = x_mat.as_flex_double_matrix()
write_2d(x)
