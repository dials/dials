def raw_2d_integration(data2d, mask2d):
    from dials.algorithms.background.flat_subtraction import flat_background_subtraction_2d
    import numpy
    n_col = numpy.size(data2d[0:1, :])
    n_row = numpy.size(data2d[:, 0:1])
    bkgr = flat_background_subtraction_2d(data2d, mask2d)
    #print 'bkgr =', bkgr
    #print data2d
    tot = 0
    for col in range(n_col):
        for row in range(n_row):
            tot += data2d[row, col]
    #print "tot =", tot
    return tot
