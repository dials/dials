import sys
import numpy
from iotbx.detectors import ImageFactory

image = ImageFactory( sys.argv[1] )
image.read()

nfast = image.parameters['SIZE1']
nslow = image.parameters['SIZE2']

print "nslow =", nslow
xsize = nslow
print "nfast =", nfast
ysize = nfast

data = image.get_raw_data()

# data2d = numpy.reshape( numpy.array( data, dtype = 'i2' ), ( ysize, xsize ) )
data2d = numpy.reshape( numpy.array( data, dtype = 'intc' ), ( ysize, xsize ) )

print "max(data2d) =", numpy.max( data2d )
print "min(data2d) =", numpy.min( data2d )

f = open( 'img_tst.raw', 'wb' )
f.write( data2d )
f.close()

'''f2 = open( 'img_raw_tst01_a.raw', 'wb' )
f2.write( numpy.transpose( data2d ) )
f2.close()'''
