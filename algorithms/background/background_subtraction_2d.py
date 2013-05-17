def curved_background_subtraction_2d(data2d, mask2d):

    import numpy
    n_col = numpy.size(data2d[0:1, :])
    n_row = numpy.size(data2d[:, 0:1])
    avg_bkgr_tot = 0.0
    avg_bkgr_cont = 0
    for row in range(0, n_row, 1):
        for col in range(0, n_col, 1):
            if mask2d[row, col] == 1:

                loc_bkgr_tot = 0.0
                loc_bkgr_cont = 0
                if mask2d[n_row - 1, col] == 0:
                    loc_bkgr_tot += data2d[n_row - 1, col]
                    loc_bkgr_cont += 1
                if mask2d[0, col] == 0:
                    loc_bkgr_tot += data2d[0, col]
                    loc_bkgr_cont += 1
                if mask2d[row, n_col - 1] == 0:
                    loc_bkgr_tot += data2d[row, n_col - 1]
                    loc_bkgr_cont += 1
                if mask2d[row, 0] == 0:
                    loc_bkgr_tot += data2d[row, 0]
                    loc_bkgr_cont += 1
                loc_bkgr = loc_bkgr_tot / float(loc_bkgr_cont)

                avg_bkgr_tot += loc_bkgr
                avg_bkgr_cont += 1

                if data2d[row, col] > loc_bkgr:
                    data2d[row, col] = data2d[row, col] - loc_bkgr
                else:
                    data2d[row, col] = 0
    avg_bkgr = avg_bkgr_tot / float(avg_bkgr_cont)
    for row in range(0, n_row, 1):
        for col in range(0, n_col, 1):
            if mask2d[row, col] == 0:
                data2d[row, col] = 0
    return avg_bkgr

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
