from __future__ import division
import numpy

data2d = numpy.arange( 36, dtype = float ).reshape( 6, 6 )
a = numpy.sum( data2d[2:4, 1:3] )
data2d[3, 3] = a
print data2d
cont = 0
for scan in data2d[2:4, 2:4]:
  cont = cont + 1
  print scan
  print cont
  scan = cont
print "data2d"
print data2d
print"______________________________________"

cont = 0
for xscan in range( 2, 4, 1 ):
  for yscan in range( 2, 4, 1 ):
    cont = cont + 1
    print cont
    print data2d[yscan, xscan]
    data2d[yscan, xscan] = cont
print "data2d"
print data2d
print"______________________________________"

a_2d = numpy.arange( 16, dtype = float ).reshape( 4, 4 )
b_2d = numpy.arange( 16, dtype = float ).reshape( 4, 4 )


print "a_2d"
print a_2d

# print "b_2d"
# print b_2d

sum_2d = numpy.sum( a_2d[1:4, 1:4] )
print "sum_2d"
print sum_2d

random_test='''
from scitbx.array_family import flex
w_var=flex.double(flex.grid(4,4))
flex.show(w_var)
'''
