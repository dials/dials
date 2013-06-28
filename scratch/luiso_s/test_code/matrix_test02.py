from __future__ import division
import numpy
import matrix_test02_sub
ysize = xsize = 16
mat2d = numpy.arange( ysize * xsize, dtype = 'uintc' ).reshape( ysize, xsize )
print  mat2d
print 'id(mat2d) =', id( mat2d )
x = 5
matrix_test02_sub.mpart( mat2d[4:8, 2:6] )
print  mat2d
print 'id(mat2d) =', id( mat2d )
