from __future__ import division
import numpy
import time
ysize = xsize = 10



mat2d1 = numpy.arange(ysize * xsize, dtype = 'uintc').reshape(ysize, xsize)
mat2d2 = numpy.transpose(mat2d1)
mat2dres = numpy.arange(ysize * xsize, dtype = 'uintc').reshape(ysize, xsize)
mat2dsmoth = numpy.zeros(ysize * xsize, dtype = 'uintc').reshape(ysize, xsize)

print time.time()
time1 = time.time()
print "time1 =", time1

for y in range(ysize):
    for x in range(xsize):
        mat2dres[y, x] = mat2d1[y, x] + mat2d2[y, x]

time2 = time.time()
print "time2 =", time2
timedif = time2 - time1
print "timedif =", timedif

#for y in range(1, ysize - 1):
#    for x in range(1, xsize - 1):
#        suma = mat2dres[y - 1, x - 1] + mat2dres[y - 1, x] + mat2dres[y - 1, x + 1] \
#             + mat2dres[y, x - 1] + mat2dres[y, x + 1]             \
#             + mat2dres[y + 1, x - 1] + mat2dres[y + 1, x] + mat2dres[y + 1, x + 1]
#        mat2dsmoth[y, x] = suma / 8

#        pscan = numpy.sum(mat2dres[y - 1:y + 2, x - 1:x + 2]) / 9.0
#        mat2dsmoth[y, x] = int(pscan)

data_tmp = numpy.copy(mat2dres[0:5, 0:5])
mat2dsmoth[0:5, 0:5] = data_tmp

data_tmp = numpy.copy(mat2dres[0:5, 5:10])
mat2dsmoth[0:5, 5:10] = data_tmp

data_tmp = numpy.copy(mat2dres[5:10, 0:5])
mat2dsmoth[5:10, 0:5] = data_tmp

data_tmp = numpy.copy(mat2dres[5:10, 5:10])
mat2dsmoth[5:10, 5:10] = data_tmp

time3 = time.time()
print "time3 =", time3
timedif = time3 - time2
print "timedif =", timedif

print 'mat2d1'
print mat2d1
print 'mat2d2'
print mat2d2
print 'mat2dres'
print mat2dres
print 'mat2dsmoth'
print mat2dsmoth
