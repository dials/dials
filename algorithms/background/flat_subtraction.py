def flat_background_subtraction_2d(data2d, diffdata2d_ext):

    import numpy
    n_col = numpy.size(data2d[0:1, :])
    n_row = numpy.size(data2d[:, 0:1])

    tot_bkgr = 0.0
    cont = 0.0
    for row in range(0, n_row, 1):
        for col in range(0, n_col, 1):
            if diffdata2d_ext[row, col] == 0:
                cont += 1
                tot_bkgr += data2d[row, col]
    if tot_bkgr > 0 and cont > 0:
        bkgr = tot_bkgr / cont
    else:
        bkgr = 0
    #print 'bkgr=', bkgr
    for row in range(0, n_row, 1):
        for col in range(0, n_col, 1):
            if diffdata2d_ext[row, col] == 1 and data2d[row, col] > bkgr:
                data2d[row, col] = data2d[row, col] - bkgr
            else:
                data2d[row, col] = 0
    return bkgr
