from __future__ import division
import sys
import numpy
from iotbx.detectors import ImageFactory
from matplotlib import pyplot as plt

image = ImageFactory( sys.argv[1] )
image.read()

nfast = image.parameters['SIZE1']
nslow = image.parameters['SIZE2']

data = image.get_raw_data()
print 'here 1'
data2d = numpy.reshape( numpy.array( data, dtype = float ), ( nfast , nslow ) )
print 'here 2'
data2dsmoth = numpy.zeros( nfast * nslow, dtype = float ).reshape( nfast, nslow )
diffdata2d = numpy.zeros( nfast * nslow, dtype = float ).reshape( nfast, nslow )

print "nslow, nfast =", nslow, nfast

print "max(data2d) =", numpy.max( data2d )
print "min(data2d) =", numpy.min( data2d )

for f in range( 1, nfast - 1 ):
    for s in range( 1, nslow - 1 ):
        pscan = float( numpy.sum( data2d[f - 1:f + 1, s - 1:s + 1] ) / 9.0 )
        data2dsmoth[f, s] = pscan

print "max(data2dsmoth) =", numpy.max( data2dsmoth )
print "min(data2dsmoth) =", numpy.min( data2dsmoth )

for f in range( 1, nfast - 1 ):
    for s in range( 1, nslow - 1 ):
        if data2d[f, s] > data2dsmoth[f, s]:
            diffdata2d[f, s] = 1

print "max(diffdata2d) =", numpy.max( diffdata2d )
print "min(diffdata2d) =", numpy.min( diffdata2d )

print "__________________________________________"

print "Plotting data2d"
plt.imshow( numpy.transpose( data2d ), interpolation = "nearest" )
plt.show()

print "Plotting data2dsmoth"
plt.imshow( numpy.transpose( data2dsmoth ), interpolation = "nearest" )
plt.show()

print "Plotting diffdata2d"
plt.imshow( numpy.transpose( diffdata2d ), interpolation = "nearest" )
plt.show()
