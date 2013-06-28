from __future__ import division
# -*- coding: utf-8 -*-
import numpy
from matplotlib import pyplot as plt
def start( data2d, xcoord, ycoord , cntrd_xcoord, cntrd_ycoord, x_sigma, y_sigma ):
    data2d = numpy.transpose( data2d )
    n_x = numpy.size( data2d[:, 0:1] )
    n_y = numpy.size( data2d[0:1, :] )
    print 'n_x =', n_x
    print 'n_y =', n_y
    data2dtmp = data2d

    data2dsmoth = numpy.zeros( n_x * n_y, dtype = float ).reshape( n_x, n_y )
    diffdata2d = numpy.zeros( n_x * n_y, dtype = int ).reshape( n_x, n_y )
    diffdata2d_ext = numpy.zeros( n_x * n_y, dtype = int ).reshape( n_x, n_y )
    paintmask = numpy.zeros( n_x * n_y, dtype = int ).reshape( n_x, n_y )
    for times in range( 5 ):
        for x in range( 1, n_x - 1 ):
            for y in range( 1, n_y - 1 ):
                pscan = float( numpy.sum( data2dtmp[x - 1:x + 2, y - 1:y + 2] ) / 9.0 )
                data2dsmoth[x, y] = pscan
        data2dtmp = data2dsmoth

    data2dsmoth = data2dsmoth + 15.0
    print "max(data2dsmoth) =", numpy.max( data2dsmoth )
    print "min(data2dsmoth) =", numpy.min( data2dsmoth )

    for x in range( 1, n_x - 1 ):
        for y in range( 1, n_y - 1 ):
            if data2d[x, y] > data2dsmoth[x, y]:
                diffdata2d[x, y] = 1

    for x in range( 1, n_x - 1 ):
        for y in range( 1, n_y - 1 ):
            if diffdata2d[x, y] == 1:
                diffdata2d_ext[x - 1:x + 2, y - 1:y + 2] = 1


    for pos in range( len( xcoord ) ):
        x_num_sum = 0.0
        y_num_sum = 0.0
        den_sum = 0.0
        for x_scan in range( int( xcoord[pos] ) - 5, int( xcoord[pos] + 6 ), 1 ):
            for y_scan in range( int( ycoord[pos] ) - 5, int( ycoord[pos] ) + 6, 1 ):
                if diffdata2d_ext[x_scan, y_scan] == 1:
                    x_num_sum = x_num_sum + float( data2d[x_scan, y_scan] ) * float( x_scan )
                    y_num_sum = y_num_sum + float( data2d[x_scan, y_scan] ) * float( y_scan )
                    den_sum = den_sum + float( data2d[x_scan, y_scan] )
        if den_sum > 0.0:
            cntrd_xcoord[pos] = x_num_sum / den_sum
            cntrd_ycoord[pos] = y_num_sum / den_sum
        else:
            print 'den_sum =', den_sum
            cntrd_xcoord[pos] = -1
            cntrd_ycoord[pos] = -1

    for pos in range( len( xcoord ) ):
        x_num_sum = 0.0
        y_num_sum = 0.0
        den_sum = 0.0
        for x_scan in range( int( xcoord[pos] ) - 5, int( xcoord[pos] ) + 6, 1 ):
            for y_scan in range( int( ycoord[pos] ) - 5, int( ycoord[pos] ) + 6, 1 ):
                if diffdata2d_ext[x_scan, y_scan] == 1:
                    x_num_sum = x_num_sum + float( data2d[x_scan, y_scan] ) * ( float( x_scan ) - float( cntrd_xcoord[pos] ) ) ** 2.0
                    y_num_sum = y_num_sum + float( data2d[x_scan, y_scan] ) * ( float( y_scan ) - float( cntrd_ycoord[pos] ) ) ** 2.0
                    den_sum = den_sum + float( data2d[x_scan, y_scan] )
        if den_sum > 0.0:
            x_sigma[pos] = numpy.sqrt( x_num_sum / den_sum )
            y_sigma[pos] = numpy.sqrt( y_num_sum / den_sum )
        else:
            print 'den_sum =', den_sum
            x_sigma[pos] = -1
            y_sigma[pos] = -1

#    print "Plotting data2dsmoth"
#    plt.imshow( numpy.transpose( data2dsmoth ), interpolation = "nearest", origin = 'lower' )
#    plt.show()

    paintmask = diffdata2d_ext
    for pos in range( len( xcoord ) ):
        paintmask[xcoord[pos] - 5:xcoord[pos] + 6, ycoord[pos] - 5:ycoord[pos] + 6] = paintmask[xcoord[pos] - 5:xcoord[pos] + 6, ycoord[pos] - 5:ycoord[pos] + 6] + 1
        nint_xcoord = int( xcoord[pos] + .5 )
        nint_ycoord = int( ycoord[pos] + .5 )
        paintmask[nint_xcoord, nint_ycoord ] = paintmask[nint_xcoord, nint_ycoord] + 1
    for pos in range( len( xcoord ) ):
        paintmask[int( cntrd_xcoord[pos] + 0.5 ), int( cntrd_ycoord[pos] + 0.5 )] = 5
    print "Plotting paintmask"
    plt.imshow( numpy.transpose( paintmask ), interpolation = "nearest", origin = 'lower' )
    plt.show()

    print "Plotting diffdata2d"
    plt.imshow( numpy.transpose( diffdata2d ), interpolation = "nearest", origin = 'lower' )
    plt.show()
#
#    print "Plotting data2d"
#    plt.imshow( numpy.transpose( data2d ), interpolation = "nearest", origin = 'lower' )
#    plt.show()
#
