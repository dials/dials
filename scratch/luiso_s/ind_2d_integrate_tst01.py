# -*- coding: utf-8 -*-
import numpy

def start( data2d, xcoord, ycoord , x_cm, y_cm, x_sigma, y_sigma ):

    data2d = numpy.transpose( data2d )

    ar_size = 8

    for pos in range( len( xcoord ) ):

        from_x = int( xcoord[pos] ) - ar_size
        x_to = int( xcoord[pos] + ar_size + 1 )
        from_y = int( ycoord[pos] - ar_size )
        y_to = int( ycoord[pos] + ar_size + 1 )
        x_cm[pos], y_cm[pos], x_sigma[pos] , y_sigma[pos] = single_spot_integrate( data2d[from_x:x_to, from_y:y_to] )
        x_cm[pos] = float( from_x ) + x_cm[pos]
        y_cm[pos] = float( from_y ) + y_cm[pos]

def single_spot_integrate( data2d ):
    n_x = numpy.size( data2d[:, 0:1] )
    n_y = numpy.size( data2d[0:1, :] )
    data2dsmoth = numpy.zeros( n_x * n_y, dtype = float ).reshape( n_x, n_y )
    diffdata2d = numpy.zeros( n_x * n_y, dtype = int ).reshape( n_x, n_y )
    diffdata2d_ext = numpy.zeros( n_x * n_y, dtype = int ).reshape( n_x, n_y )
    data2dtmp = data2d

    x_to = n_x
    y_to = n_y

    threshold_shift = 7.0
    ext_area = 4

    for times in range( 5 ):
        for x in range( 0, x_to, 1 ):
            for y in range( 0, y_to, 1 ):
                pscan = float( numpy.sum( data2dtmp[x - 1:x + 2, y - 1:y + 2] ) / 9.0 )
                data2dsmoth[x, y] = pscan
    data2dtmp = data2dsmoth

    data2dsmoth[0:x_to, 0:y_to] = data2dsmoth[0:x_to, 0:y_to] + threshold_shift

    for x in range( 0, x_to, 1 ):
        for y in range( 0, y_to, 1 ):
            if data2d[x, y] > data2dsmoth[x, y]:
                diffdata2d[x, y] = 1

    for x in range( 0, x_to, 1 ):
        for y in range( 0, y_to, 1 ):
            if diffdata2d[x, y] == 1:
                diffdata2d_ext[x - ext_area:x + ext_area + 1, y - ext_area:y + ext_area + 1] = 1

    x_num_sum = 0.0
    y_num_sum = 0.0
    den_sum = 0.0
    for x in range( 0, x_to, 1 ):
        for y in range( 0, y_to, 1 ):
            if diffdata2d_ext[x, y] == 1:
                x_num_sum = x_num_sum + float( data2d[x, y] ) * float( x )
                y_num_sum = y_num_sum + float( data2d[x, y] ) * float( y )
                den_sum = den_sum + float( data2d[x, y] )
    if den_sum > 0.0:
        x_cm = x_num_sum / den_sum
        y_cm = y_num_sum / den_sum
    else:
        print 'den_sum =', den_sum
        x_cm = -1
        y_cm = -1

    x_num_sum = 0.0
    y_num_sum = 0.0
    den_sum = 0.0
    for x in range( 0, x_to, 1 ):
        for y in range( 0, y_to, 1 ):
            if diffdata2d_ext[x, y] == 1:
                x_num_sum = x_num_sum + float( data2d[x, y] ) * ( float( x ) - float( x_cm ) ) ** 2.0
                y_num_sum = y_num_sum + float( data2d[x, y] ) * ( float( y ) - float( y_cm ) ) ** 2.0
                den_sum = den_sum + float( data2d[x, y] )
    if den_sum > 0.0:
        x_sg = numpy.sqrt( x_num_sum / den_sum )
        y_sg = numpy.sqrt( y_num_sum / den_sum )
    else:
        print 'den_sum =', den_sum
        x_sg = -1
        y_sg = -1

    return x_cm, y_cm, x_sg, y_sg
