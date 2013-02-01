import sys
import numpy
from iotbx.detectors import ImageFactory
from matplotlib import pyplot as plt

image = ImageFactory(sys.argv[1])
image.read()

nfast = image.parameters['SIZE1']
nslow = image.parameters['SIZE2']

data = image.get_raw_data()
data2d = numpy.arange(nfast * nslow).reshape(nfast,nslow)
data2dsmoth = numpy.arange(nfast * nslow).reshape(nfast,nslow)
diffdata2d = numpy.arange(nfast * nslow).reshape(nfast,nslow)

for f in range(nfast):
    for s in range(nslow):
        data2d[f,s] = data[s * nfast + f]

print "nslow, nfast =", nslow, nfast

print "max(data2d) =", numpy.max(data2d)
print "min(data2d) =", numpy.min(data2d)


for f in range(1,nfast-1):
    for s in range(1,nslow-1):
        data2dsmoth[f,s]=( data2d[f-1,s] + data2d[f+1,s] + data2d[f,s-1] + data2d[f,s+1] ) / 4

print "max(data2dsmoth) =", numpy.max(data2dsmoth)
print "min(data2dsmoth) =", numpy.min(data2dsmoth)



for f in range(1,nfast-1):
    for s in range(1,nslow-1):
        if diffdata2d[f,s] < 0:
            diffdata2d[f,s] = 0

diffdata2d = data2d - data2dsmoth
print "max(diffdata2d) =", numpy.max(diffdata2d)
print "min(diffdata2d) =", numpy.min(diffdata2d)


plt.imshow(data2d)
plt.savefig('img_data2d.png')

plt.imshow(diffdata2d)
plt.savefig('img_diff.png')

plt.imshow(data2dsmoth)
plt.savefig('img_smoth.png')