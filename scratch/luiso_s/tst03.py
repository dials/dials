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
data2d = numpy.arange(ysize * xsize).reshape(ysize,xsize)
data2dsmoth = numpy.arange(ysize * xsize).reshape(ysize,xsize)
diffdata2d = numpy.arange(ysize * xsize).reshape(ysize,xsize)



for y in range(ysize):
    for x in range(xsize):
        data2d[y,x] = data[x * ysize + y]

print "max(data2d) =", numpy.max(data2d)
print "min(data2d) =", numpy.min(data2d)

for y in range(1, ysize-1):
    for x in range(1, xsize-1):
        data2dsmoth[y, x]=( data2d[y-1, x] + data2d[y+1, x] + data2d[y, x-1] + data2d[y, x+1] ) / 4

print "max(data2dsmoth) =", numpy.max(data2dsmoth)
print "min(data2dsmoth) =", numpy.min(data2dsmoth)

data2dsmoth = data2dsmoth + 5

for y in range(1,ysize-1):
    for x in range(1,xsize-1):
        if data2d[y, x] > data2dsmoth[y, x]:
            diffdata2d[y, x] = 100


diffdata2d = 0

print "max(diffdata2d) =", numpy.max(diffdata2d)
print "min(diffdata2d) =", numpy.min(diffdata2d)


plt.imshow(data2d)
plt.savefig('img_data2d.png')

plt.imshow(diffdata2d)
plt.savefig('img_diff.png')

plt.imshow(data2dsmoth)
plt.savefig('img_smoth.png')

