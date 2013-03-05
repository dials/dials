import numpy
from matplotlib import pyplot as plt
def start( data2d, xcoord, ycoord , cntrd_xcoord, cntrd_ycoord, pos_sigma ):
    n_x = numpy.size( data2d[:, 0:1] )
    n_y = numpy.size( data2d[0:1, :] )
    print 'n_x =', n_x
    print 'n_y =', n_y
    data2dtmp = data2d

    data2dsmoth = numpy.zeros( n_x * n_y, dtype = float ).reshape( n_x, n_y )
    diffdata2d = numpy.zeros( n_x * n_y, dtype = float ).reshape( n_x, n_y )
    paintmask = numpy.zeros( n_x * n_y, dtype = float ).reshape( n_x, n_y )



    for times in range( 5 ):
        for x in range( 3, n_x - 3 ):
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
                paintmask[x, y] = 1

    print "max(diffdata2d) =", numpy.max( diffdata2d )
    print "min(diffdata2d) =", numpy.min( diffdata2d )
    print "__________________________________________"

    # for coord in coordlist:
    for pos in range( len( xcoord ) ):
        paintmask[xcoord[pos] - 5:xcoord[pos] + 6, ycoord[pos] - 5:ycoord[pos] + 6] = paintmask[xcoord[pos] - 5:xcoord[pos] + 6, ycoord[pos] - 5:ycoord[pos] + 6] + 1
        paintmask[xcoord[pos], ycoord[pos]] = paintmask[xcoord[pos], ycoord[pos]] + 1

    for pos in range( len( xcoord ) ):
        x_num_sum = 0.0
        y_num_sum = 0.0
        den_sum = 0.0
        for x_scan in range( xcoord[pos] - 5, xcoord[pos] + 6, 1 ):
            for y_scan in range( ycoord[pos] - 5, ycoord[pos] + 6, 1 ):
                if diffdata2d[x_scan, y_scan] == 1:
                    x_num_sum = x_num_sum + data2d[x_scan, y_scan] * x_scan
                    y_num_sum = y_num_sum + data2d[x_scan, y_scan] * y_scan
                    den_sum = den_sum + data2d[x_scan, y_scan]
        cntrd_xcoord[pos] = x_num_sum / den_sum
        cntrd_ycoord[pos] = y_num_sum / den_sum
        paintmask[int( cntrd_xcoord[pos] + .5 ), int( cntrd_ycoord[pos] + .5 )] = 5


    for pos in range( len( xcoord ) ):
        num_sum = 0.0
        den_sum = 0.0
        for x_scan in range( xcoord[pos] - 5, xcoord[pos] + 6, 1 ):
            for y_scan in range( ycoord[pos] - 5, ycoord[pos] + 6, 1 ):
                if diffdata2d[x_scan, y_scan] == 1:
                    num_sum = num_sum + data2d[x_scan, y_scan] * ( x_scan - cntrd_xcoord[pos] ) ** 2.0
                    den_sum = den_sum + data2d[x_scan, y_scan]
        pos_sigma[pos] = numpy.sqrt( num_sum / den_sum )


    print "Plotting data2dsmoth"
    plt.imshow( numpy.transpose( data2dsmoth ), interpolation = "nearest" )
    plt.show()

    print "Plotting paintmask"
    plt.imshow( numpy.transpose( paintmask ), interpolation = "nearest" )
    plt.show()

    print "Plotting diffdata2d"
    plt.imshow( numpy.transpose( diffdata2d ), interpolation = "nearest" )
    plt.show()

    print "Plotting data2d"
    plt.imshow( numpy.transpose( data2d ), interpolation = "nearest" )
    plt.show()
