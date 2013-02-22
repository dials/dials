import numpy
import time
ysize = xsize = 2000



mat2d1 = numpy.arange( ysize * xsize, dtype = 'uintc' ).reshape( ysize, xsize )
mat2d2 = numpy.transpose( mat2d1 )
mat2dres = numpy.arange( ysize * xsize, dtype = 'uintc' ).reshape( ysize, xsize )

print time.time()
time1 = time.time()
print "time1 =", time1

for y in range( 1, ysize - 1 ):
    for x in range( 1, xsize - 1 ):
        mat2dres[y, x] = mat2d1[y, x] + mat2d2[y, x]

time2 = time.time()
print "time2 =", time2
timedif = time2 - time1
print "timedif =", timedif
