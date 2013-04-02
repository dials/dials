import numpy

def funt():
    print 'hola por aqui tst 01'



def find_mask_2d(data2d):

    n_col = numpy.size(data2d[0:1, :])
    n_row = numpy.size(data2d[:, 0:1])
    #print 'n_col,n_row =', n_col, n_row

    data2dsmoth = numpy.zeros(n_row * n_col, dtype = float).reshape(n_row, n_col)

    diffdata2d = numpy.zeros(n_row * n_col, dtype = int).reshape(n_row, n_col)
    diffdata2d_ext = numpy.zeros(n_row * n_col, dtype = int).reshape(n_row, n_col)
    data2dtmp = data2d

    for times in range(5):
        for row in range(1, n_row - 1, 1):
            for col in range(1, n_col - 1, 1):
                pscan = float(numpy.sum(data2dtmp[row - 1:row + 2, col - 1:col + 2]) / 9.0)
                data2dsmoth[row, col] = pscan
    data2dtmp = data2dsmoth
    #threshold_shift = (numpy.max(data2dsmoth) - numpy.min(data2dsmoth)) / 2.0 # This used to be one of this "magical variables"

#######################################################################################################
    cont = 0                                                                  # This way to calculate
    dif_tot = 0                                                               # this magical variable
    for row in range(0, n_row, 1):                                            # is more statistical
        for col in range(0, n_col, 1):                                        # and seems to be giving
            cont += 1                                                         # better results
            dif_tot += numpy.abs(data2d[row, col] - data2dsmoth[row, col])    #
    dif_avg = dif_tot / cont                                                  #
    #print 'dif_avg=', dif_avg                                                #
    threshold_shift = dif_avg * 3.0                                           #
#######################################################################################################
    print 'threshold_shift =', threshold_shift

    data2dsmoth[0:n_row, 0:n_col] = data2dsmoth[0:n_row, 0:n_col] + threshold_shift
    ext_area = 1                                                               # This used to be one of this "magical variables"
    for row in range(0, n_row, 1):
        for col in range(0, n_col, 1):
            if data2d[row, col] > data2dsmoth[row, col]:
                diffdata2d[row, col] = 1

    for row in range(ext_area, n_row - ext_area + 1, 1):
        for col in range(ext_area, n_col - ext_area + 1, 1):
            if diffdata2d[row, col] == 1:
                diffdata2d_ext[row - ext_area:row + ext_area + 1, col - ext_area:col + ext_area + 1] = 1

    #print '_____________diffdata2d'
    #print diffdata2d
    #print '_____________diffdata2d_ext'
    #print diffdata2d_ext
    return diffdata2d_ext
'''
############################################################################### flat background
    tot_bkgr = 0.0                                                            # version
    cont = 0.0                                                                #
    for row in range(0, n_row, 1):                                            #
        for col in range(0, n_col, 1):                                        #
            if diffdata2d_ext[row, col] == 0:                                 #
                cont += 1                                                     #
                tot_bkgr += data2d[row, col]                                  #
    bkgr = tot_bkgr / cont                                                    #
    #print 'bkgr=', bkgr                                                      #
    for row in range(0, n_row, 1):                                            #
        for col in range(0, n_col, 1):                                        #
            if diffdata2d_ext[row, col] == 1 and data2d[row, col] > bkgr:     #
                data2d[row, col] = data2d[row, col] - bkgr                    #
            else:                                                             #
                data2d[row, col] = 0                                          #
###############################################################################
'''
############################################################################### curved background
#   colbord = int(n_col / 5)                                                     # version
#   colbord = int(n_row / 5)                                                     #
#   for row in range(0, n_row, 1):                                               #
#       for col in range(0, n_col, 1):                                           #
#           if diffdata2d_ext[row,col] == 1:                                     #
#               top_av = float(numpy.sum(data2d[:colbord, col - 1:col + 2]))        #
#               bot_av = float(numpy.sum(data2d[n_row - colbord:, col - 1:col + 2])) #
#               lft_av = float(numpy.sum(data2d[row - 1:row + 2, colbord]))         #
#               rgt_av = float(numpy.sum(data2d[row - 1:row + 2, n_col - colbord:])) #
#               bkgr = (top_av + bot_av + lft_av + rgt_av) / 4.0              #
#               if data2d[row,col] > bkgr:                                       #
#                   data2d[row,col] = data2d[row,col] - bkgr                        #
#               else:                                                         #
#                   data2d[row,col] = 0                                          #
#   for row in range(0, n_row, 1):                                               #
#       for col in range(0, n_col, 1):                                           #
#           if diffdata2d_ext[row,col] == 0:                                     #
#               data2d[row,col] = 0                                              #
###############################################################################

#   #col_num_sum = 0.0
#   #row_num_sum = 0.0
#   itst_sum = 0.0
#   for row in range(0, n_row, 1):
#       for col in range(0, n_col, 1):
#           itst_sum += float(data2d[row, col])
#   if itst_sum > 0.0:
#       tot = itst_sum
#   else:
#       print 'itst_sum =', itst_sum
#       tot = -1
#
##   display_image_with_predicted_spots_n_centoids(data2d, col_cm, row_cm, n_col / 2, n_row / 2)
#
#   return tot
