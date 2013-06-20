#
# Background Subtraction 2D lower level numpy version
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
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

                if loc_bkgr_cont > 0:
                    loc_bkgr = loc_bkgr_tot / float(loc_bkgr_cont)
                else:
                    loc_bkgr = 0
                avg_bkgr_tot += loc_bkgr
                avg_bkgr_cont += 1

                if data2d[row, col] > loc_bkgr:
                    data2d[row, col] = data2d[row, col] - loc_bkgr
                else:
                    data2d[row, col] = 0
    if avg_bkgr_cont > 0:
        avg_bkgr = avg_bkgr_tot / float(avg_bkgr_cont)
    else:
        avg_bkgr = 0
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
        avg_bkgr = tot_bkgr / cont
    else:
        avg_bkgr = 0
    #print 'avg_bkgr=', avg_bkgr
    for row in range(0, n_row, 1):
        for col in range(0, n_col, 1):
            if diffdata2d_ext[row, col] == 1 and data2d[row, col] > avg_bkgr:
                data2d[row, col] = data2d[row, col] - avg_bkgr
            else:
                data2d[row, col] = 0

    return avg_bkgr

def flat_background_calc_2d(data2d, diffdata2d_ext):

    import numpy
    n_col = numpy.size(data2d[0:1, :])
    n_row = numpy.size(data2d[:, 0:1])
    avg_bkgr_2d = numpy.copy(data2d)

    tot_bkgr = 0.0
    cont = 0.0
    for row in range(0, n_row, 1):
        for col in range(0, n_col, 1):
            if diffdata2d_ext[row, col] == 0:
                cont += 1
                tot_bkgr += data2d[row, col]
    if tot_bkgr > 0 and cont > 0:
        avg_bkgr = tot_bkgr / cont
    else:
        avg_bkgr = 0
    #print 'avg_bkgr=', avg_bkgr
    for row in range(0, n_row, 1):
        for col in range(0, n_col, 1):
            if diffdata2d_ext[row, col] != 0:
                avg_bkgr_2d[row, col] = avg_bkgr

            #if diffdata2d_ext[row, col] == 1 and data2d[row, col] > avg_bkgr:
            #    data2d[row, col] = data2d[row, col] - avg_bkgr
            #else:
            #    data2d[row, col] = 0
    return avg_bkgr_2d

def curved_background_calc_2d(data2d, mask2d):

    import numpy
    n_col = numpy.size(data2d[0:1, :])
    n_row = numpy.size(data2d[:, 0:1])
    avg_bkgr_2d = numpy.copy(data2d)
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

                if loc_bkgr_cont > 0:
                    loc_bkgr = loc_bkgr_tot / float(loc_bkgr_cont)
                else:
                    loc_bkgr = 0

                avg_bkgr_2d[row, col] = loc_bkgr


    return avg_bkgr_2d

