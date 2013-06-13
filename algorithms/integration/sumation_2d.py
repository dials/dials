def raw_2d_integration_old(data2d, mask2d, bkgr):
    import numpy
    n_col = numpy.size(data2d[0:1, :])
    n_row = numpy.size(data2d[:, 0:1])

    print 'bkgr =', bkgr
    #print data2d
    i_tot = 0
    npix_bkgr = 0
    npix_mask = 0
    for col in range(n_col):
        for row in range(n_row):

            if mask2d[row, col] == 1 :
                i_tot += data2d[row, col]
                npix_mask += 1
            else:
                npix_bkgr += 1
    sig = numpy.sqrt(i_tot + (1.0 + (npix_mask) / (npix_bkgr)) * (npix_mask * bkgr))
    #print "i_tot =", i_tot
    return i_tot, sig


def raw_2d_integration(data2d, mask2d, bkgr2d):
    import numpy
    n_col = numpy.size(data2d[0:1, :])
    n_row = numpy.size(data2d[:, 0:1])

    #print data2d
    i_tot = 0
    npix_bkgr = 0
    npix_mask = 0
    cont = 0
    tot_bkgr = 0
    for col in range(n_col):
        for row in range(n_row):

            if mask2d[row, col] == 1 :
                i_tot = i_tot + data2d[row, col] - bkgr2d[row, col]
                npix_mask += 1
            else:
                npix_bkgr += 1

            cont += 1
            tot_bkgr += bkgr2d[row, col]
    if tot_bkgr > 0 and cont > 0:
        bkgr = tot_bkgr / cont
    else:
        bkgr = 0
    print 'bkgr =', bkgr

    sig = numpy.sqrt(i_tot + (1.0 + (npix_mask) / (npix_bkgr)) * (npix_mask * bkgr))
    #print "i_tot =", i_tot
    return i_tot, sig
