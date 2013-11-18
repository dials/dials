from __future__ import division
import sys
import numpy
from iotbx.detectors import ImageFactory
from matplotlib import pyplot as plt

image = ImageFactory(sys.argv[1])
image.read()

nfast = image.parameters['SIZE1']
nslow = image.parameters['SIZE2']

print "nslow =", nslow
xsize = nslow
print "nfast =", nfast
ysize = nfast

data = image.get_raw_data()
data2d = numpy.reshape(numpy.array(data, dtype = 'uintc'),(ysize, xsize))
data2dsmoth = numpy.zeros(nfast * nslow, dtype='uintc').reshape(ysize, xsize)
diffdata2d = numpy.zeros(nfast * nslow, dtype='uintc').reshape(ysize, xsize)

for y in range(1, ysize-1):
  for x in range(1, xsize-1):
    data2dsmoth[y, x]=( data2d[y-1, x] + data2d[y+1, x] + data2d[y, x-1] + data2d[y, x+1] ) / 4

#for y in range(1,ysize-1):
#    for x in range(1,xsize-1):
#        if data2d[y, x] > data2dsmoth[y, x]:
#            diffdata2d[y, x] = 100
#        else:
#            diffdata2d[y, x] = 10

plt.imshow(data2d)
plt.show()

plt.imshow(data2dsmoth)
plt.show()

plt.imshow(diffdata2d)
plt.show()
