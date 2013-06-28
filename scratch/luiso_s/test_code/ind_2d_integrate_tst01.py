from __future__ import division
# -*- coding: utf-8 -*-
import numpy
#path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))      # start copy-pasted from example
#if not path in sys.path:                                                   #
#    sys.path.insert(1, path)                                               #
#del path                                                                   # end copy-pasted from example
def start(data2d, xcoord, ycoord , x_cm, y_cm, x_sigma, y_sigma):

    ar_size = 8

    for pos in range(len(xcoord)):

        from_x = int(xcoord[pos]) - ar_size
        x_to = int(xcoord[pos] + ar_size + 1)
        from_y = int(ycoord[pos] - ar_size)
        y_to = int(ycoord[pos] + ar_size + 1)

        # Next line should be used as template when you want to calculate centoids and sigmas from a 2D portion of image

        x_cm[pos], y_cm[pos], x_sigma[pos] , y_sigma[pos], tot_itst = single_spot_integrate_2d(data2d[from_y:y_to , from_x:x_to])

        x_cm[pos] = float(from_x) + x_cm[pos]
        y_cm[pos] = float(from_y) + y_cm[pos]

def single_spot_integrate_2d(data2d):
    x_to = numpy.size(data2d[0:1, :])
    y_to = numpy.size(data2d[:, 0:1])
    data2dsmoth = numpy.zeros(y_to * x_to, dtype = float).reshape(y_to, x_to)

    diffdata2d = numpy.zeros(y_to * x_to, dtype = int).reshape(y_to, x_to)
    diffdata2d_ext = numpy.zeros(y_to * x_to, dtype = int).reshape(y_to, x_to)

    data2dtmp = data2d

    threshold_shift = 7.0
    ext_area = 5

    for times in range(5):
        for y in range(0, y_to, 1):
            for x in range(0, x_to, 1):
                pscan = float(numpy.sum(data2dtmp[y - 1:y + 2, x - 1:x + 2]) / 9.0)
                data2dsmoth[y, x] = pscan
        data2dtmp = data2dsmoth

    data2dsmoth[0:y_to, 0:x_to] = data2dsmoth[0:y_to, 0:x_to] + threshold_shift

    for y in range(0, y_to, 1):
        for x in range(0, x_to, 1):
            if data2d[y, x] > data2dsmoth[y, x]:
                diffdata2d[y, x] = 1

    for y in range(0, y_to, 1):
        for x in range(0, x_to, 1):
            if diffdata2d[y, x] == 1:
                diffdata2d_ext[y - ext_area:y + ext_area + 1, x - ext_area:x + ext_area + 1] = 1

    x_num_sum = 0.0
    y_num_sum = 0.0
    den_sum = 0.0
    for y in range(0, y_to, 1):
        for x in range(0, x_to, 1):
            if diffdata2d_ext[y, x] == 1:
                x_num_sum = x_num_sum + float(data2d[y, x]) * float(x)
                y_num_sum = y_num_sum + float(data2d[y, x]) * float(y)
                den_sum = den_sum + float(data2d[y, x])
    if den_sum > 0.0:
        col_cm = x_num_sum / den_sum
        row_cm = y_num_sum / den_sum
    else:
        print 'den_sum =', den_sum
        col_cm = -1
        row_cm = -1

    x_num_sum = 0.0
    y_num_sum = 0.0
#    den_sum = 0.0
    for y in range(0, y_to, 1):
        for x in range(0, x_to, 1):
            if diffdata2d_ext[y, x] == 1:
                x_num_sum = x_num_sum + float(data2d[y, x]) * (float(x) - float(col_cm)) ** 2.0
                y_num_sum = y_num_sum + float(data2d[y, x]) * (float(y) - float(row_cm)) ** 2.0
#                den_sum = den_sum + float(data2d[y, x])
    if den_sum > 0.0:
        col_sig = numpy.sqrt(x_num_sum / den_sum)
        row_sig = numpy.sqrt(y_num_sum / den_sum)
    else:
        print 'den_sum =', den_sum
        col_sig = -1
        row_sig = -1
    return row_cm, col_cm, row_sig, col_sig, den_sum



#def single_spot_integrate( data2d ):
#    x_to = numpy.size( data2d[0:1, :] )
#    y_to = numpy.size( data2d[:, 0:1] )
#    data2dsmoth = numpy.zeros( y_to * x_to, dtype = float ).reshape( y_to, x_to )
#
#    diffdata2d = numpy.zeros( y_to * x_to, dtype = int ).reshape( y_to, x_to )
#    diffdata2d_ext = numpy.zeros( y_to * x_to, dtype = int ).reshape( y_to, x_to )
#
#    data2dtmp = data2d
#
#    threshold_shift = 7.0
#    ext_area = 4
#
#    for times in range( 5 ):
#        for y in range( 0, y_to, 1 ):
#            for x in range( 0, x_to, 1 ):
#                pscan = float( numpy.sum( data2dtmp[y - 1:y + 2, x - 1:x + 2] ) / 9.0 )
#                data2dsmoth[y, x] = pscan
#        data2dtmp = data2dsmoth
#
#    data2dsmoth[0:y_to, 0:x_to] = data2dsmoth[0:y_to, 0:x_to] + threshold_shift
#
#    for y in range( 0, y_to, 1 ):
#        for x in range( 0, x_to, 1 ):
#            if data2d[y, x] > data2dsmoth[y, x]:
#                diffdata2d[y, x] = 1
#
#    for y in range( 0, y_to, 1 ):
#        for x in range( 0, x_to, 1 ):
#            if diffdata2d[y, x] == 1:
#                diffdata2d_ext[y - ext_area:y + ext_area + 1, x - ext_area:x + ext_area + 1] = 1
#
#    x_num_sum = 0.0
#    y_num_sum = 0.0
#    den_sum = 0.0
#    for y in range( 0, y_to, 1 ):
#        for x in range( 0, x_to, 1 ):
#            if diffdata2d_ext[y, x] == 1:
#                x_num_sum = x_num_sum + float( data2d[y, x] ) * float( x )
#                y_num_sum = y_num_sum + float( data2d[y, x] ) * float( y )
#                den_sum = den_sum + float( data2d[y, x] )
#    if den_sum > 0.0:
#        x_cm = x_num_sum / den_sum
#        y_cm = y_num_sum / den_sum
#    else:
#        print 'den_sum =', den_sum
#        x_cm = -1
#        y_cm = -1
#
#    x_num_sum = 0.0
#    y_num_sum = 0.0
#    den_sum = 0.0
#    for y in range( 0, y_to, 1 ):
#        for x in range( 0, x_to, 1 ):
#            if diffdata2d_ext[y, x] == 1:
#                x_num_sum = x_num_sum + float( data2d[y, x] ) * ( float( x ) - float( x_cm ) ) ** 2.0
#                y_num_sum = y_num_sum + float( data2d[y, x] ) * ( float( y ) - float( y_cm ) ) ** 2.0
#                den_sum = den_sum + float( data2d[y, x] )
#    if den_sum > 0.0:
#        x_sg = numpy.sqrt( x_num_sum / den_sum )
#        y_sg = numpy.sqrt( y_num_sum / den_sum )
#    else:
#        print 'den_sum =', den_sum
#        x_sg = -1
#        y_sg = -1
#
#    return x_cm, y_cm, x_sg, y_sg
#
