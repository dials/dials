import numpy
import integrate_2d

from iotbx.detectors import ImageFactory
# from matplotlib import pyplot as plt

image = ImageFactory( 'thaumatin_die_M1S5_1_0001.img' )
image.read()

nfast = image.parameters['SIZE1']
nslow = image.parameters['SIZE2']

data = image.get_raw_data()
print '01 nfast, nslow =', nfast, nslow
data2d = numpy.reshape( numpy.array( data, dtype = float ), ( nfast , nslow ) )

f_from = int( nfast * 2 / 6 )  #   temporarily using
f_to = int( nfast / 2 )  #         only part of the image
s_from = int( nslow * 2 / 6 )  #   so I can code a bit faster
s_to = int( nslow / 2 )  #
nfast = f_to - f_from  #           It is taking way to long if we
nslow = s_to - s_from  #           use the entire image
data2d = data2d[f_from:f_to, s_from:f_to]
print 'new nfast, nslow =', nfast, nslow

# coordlist = [[128, 217], [139, 227]  , [147, 235], [157, 245], [166, 253]]

xcoord = [128, 139, 147, 157, 166]
ycoord = [217, 227, 235, 245, 253]

pos_sigma = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0]
cntrd_xcoord = [-1.0, -1.0, -1.0, -1.0, -1.0]
cntrd_ycoord = [-1.0, -1.0, -1.0, -1.0, -1.0]

print "nslow, nfast =", nslow, nfast
print "max(data2d) =", numpy.max( data2d )
print "min(data2d) =", numpy.min( data2d )

integrate_2d.start( data2d, xcoord, ycoord , cntrd_xcoord, cntrd_ycoord, pos_sigma )

for i in range( len( xcoord ) ):
    print 'x,y (centroid) =', cntrd_xcoord[i], cntrd_ycoord[i]
    print 'sigma (pos) =', pos_sigma [i]



