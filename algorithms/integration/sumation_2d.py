def raw_2d_integration(data2d, mask2d):
    import numpy
    n_col = numpy.size(data2d[0:1, :])
    n_row = numpy.size(data2d[:, 0:1])

    #print 'bkgr =', bkgr
    #print data2d
    tot = 0
    for col in range(n_col):
        for row in range(n_row):
            tot += data2d[row, col]
    sig = numpy.sqrt(tot)
    #print "tot =", tot
    return tot, sig
