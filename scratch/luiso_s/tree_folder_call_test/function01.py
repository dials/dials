import numpy

def funt():
    print 'hola por aqui tst 01'


def find_mask_2d(data2d):

    n_col = numpy.size(data2d[0:1, :])
    n_row = numpy.size(data2d[:, 0:1])
    #print 'n_col,n_row =', n_col, n_row
    #
    #data2dsmoth = numpy.zeros_like(data2d)
    #diffdata2d = numpy.zeros_like(data2d)
    #data2dtmp = numpy.empty_like(data2d)
    #data2dtmp[:, :] = data2d[:, :]

    data2dsmoth = numpy.zeros_like(data2d)
    diffdata2d = numpy.zeros_like(data2d)
    data2dtmp = numpy.copy(data2d)

    for times in range(5):
        for row in range(1, n_row - 1, 1):
            for col in range(1, n_col - 1, 1):
                pscan = numpy.sum(data2dtmp[row - 1:row + 2, col - 1:col + 2]) / 9.0
                data2dsmoth[row, col] = int(pscan)
        data2dtmp[:, :] = data2dsmoth[:, :]


#######################################################################################################
    #cont = 0                                                                  # This way to calculate
    #dif_tot = 0                                                               # this magical variable
    #for row in range(0, n_row, 1):                                            # is more statistical
    #    for col in range(0, n_col, 1):                                        # and seems to be giving
    #        cont += 1                                                         # better results
    #        dif_tot += numpy.abs(data2d[row, col] - data2dsmoth[row, col])    #
    #dif_avg = dif_tot / cont                                                  #
    ##print 'dif_avg=', dif_avg                                                #
    ##threshold_shift = 7.39432533334
#######################################################################################################
    threshold_shift = 7

    #print 'threshold_shift =', threshold_shift

    data2dsmoth[:, :] = data2dsmoth[:, :] + threshold_shift

    for row in range(0, n_row, 1):
        for col in range(0, n_col, 1):
            if data2d[row, col] > data2dsmoth[row, col]:
                diffdata2d[row, col] = 1

    return diffdata2d, data2dsmoth


'''
def find_mask_2d(data2d):

    n_col = numpy.size(data2d[0:1, :])
    n_row = numpy.size(data2d[:, 0:1])
    #print 'n_col,n_row =', n_col, n_row

    data2dsmoth = numpy.zeros(n_row * n_col, dtype = int).reshape(n_row, n_col)
    diffdata2d = numpy.zeros(n_row * n_col, dtype = int).reshape(n_row, n_col)
    data2dtmp = numpy.zeros(n_row * n_col, dtype = int).reshape(n_row, n_col)
    data2dtmp[:, :] = data2d[:, :]

    for times in range(5):
        for row in range(1, n_row - 1, 1):
            for col in range(1, n_col - 1, 1):
                pscan = float(numpy.sum(data2dtmp[row - 1:row + 2, col - 1:col + 2]) / 9.0)
                data2dsmoth[row, col] = int(pscan)
        data2dtmp[:, :] = data2dsmoth[:, :]


#######################################################################################################
    #cont = 0                                                                  # This way to calculate
    #dif_tot = 0                                                               # this magical variable
    #for row in range(0, n_row, 1):                                            # is more statistical
    #    for col in range(0, n_col, 1):                                        # and seems to be giving
    #        cont += 1                                                         # better results
    #        dif_tot += numpy.abs(data2d[row, col] - data2dsmoth[row, col])    #
    #dif_avg = dif_tot / cont                                                  #
    ##print 'dif_avg=', dif_avg                                                #
    ##threshold_shift = 7.39432533334
#######################################################################################################
    threshold_shift = 7

    #print 'threshold_shift =', threshold_shift

    data2dsmoth[0:n_row, 0:n_col] = data2dsmoth[0:n_row, 0:n_col] + threshold_shift

    for row in range(0, n_row, 1):
        for col in range(0, n_col, 1):
            if data2d[row, col] > data2dsmoth[row, col]:
                diffdata2d[row, col] = 1

    return diffdata2d
'''
def find_ext_mask_3d(diffdata3d):
    n_frm = numpy.size(diffdata3d[:, 0:1, 0:1])
    n_row = numpy.size(diffdata3d[0:1, :, 0:1])
    n_col = numpy.size(diffdata3d[0:1, 0:1, :])
    ext_area = 1
    diffdata3d_ext = numpy.zeros(n_row * n_col * n_frm , dtype = int).reshape(n_frm, n_row, n_col)
    for frm in range(ext_area, n_frm - ext_area + 1, 1):
        for row in range(ext_area, n_row - ext_area + 1, 1):
            for col in range(ext_area, n_col - ext_area + 1, 1):
                if diffdata3d[frm, row, col] == 1:
                    diffdata3d_ext[frm - ext_area:frm + ext_area + 1, row - ext_area:row + ext_area + 1, col - ext_area:col + ext_area + 1] = 1
    return diffdata3d_ext

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

def find_bound_2d(mask):
    n_col = numpy.size(mask[0:1, :])
    n_row = numpy.size(mask[:, 0:1])
    tmp_mask = numpy.zeros(n_row * n_col, dtype = int).reshape(n_row, n_col)
    tmp_mask[:, :] = mask[:, :]
    x_from_lst = []
    x_to_lst = []
    y_from_lst = []
    y_to_lst = []

    for row in range(0, n_row, 1):
        for col in range(0, n_col, 1):
            if mask[row, col] == 1 and tmp_mask[row, col] == 1:
                lft_bound = col - 1
                rgt_bound = col + 1
                bot_bound = row - 1
                top_bound = row + 1

                bot_old = 0
                top_old = 0
                lft_old = 0
                rgt_old = 0
                in_img = "True"
                while in_img == "True":

                    # left wall
                    if lft_bound > 0 and top_bound + 1 < n_row and bot_bound > 0:
                        for scan_row in range(bot_bound, top_bound + 1, 1):
                            if mask[scan_row, lft_bound] == 1:
                                lft_bound -= 1
                                break
                    else:
                        in_img = "false"
                        break

                    # right wall
                    if rgt_bound < n_col and top_bound + 1 < n_row and bot_bound > 0:
                        for scan_row in range(bot_bound, top_bound + 1, 1):
                            if mask[scan_row, rgt_bound] == 1:
                                rgt_bound += 1
                                break
                    else:
                        in_img = "false"
                        break
                    # bottom wall
                    if bot_bound > 0 and lft_bound > 0 and rgt_bound + 1 < n_col:
                        for scan_col in range(lft_bound, rgt_bound + 1 , 1):
                            if mask[bot_bound, scan_col] == 1:
                                bot_bound -= 1
                                break
                    else:
                        in_img = "false"
                        break
                    # top wall
                    if top_bound < n_row and lft_bound > 0 and rgt_bound + 1 < n_col:
                        for scan_col in range(lft_bound, rgt_bound + 1 , 1):
                            if mask[top_bound, scan_col] == 1:
                                top_bound += 1
                                break
                    else:
                        in_img = "false"
                        break
                    if bot_bound == bot_old and top_bound == top_old and lft_bound == lft_old and rgt_bound == rgt_old:
                        break
                    bot_old = bot_bound
                    top_old = top_bound
                    lft_old = lft_bound
                    rgt_old = rgt_bound

                #lst_coord.append([bot_bound, lft_bound, top_bound, rgt_bound])
                tmp_mask[bot_bound:top_bound, lft_bound:rgt_bound] = 0
                if  in_img == "True":
                    x_from_lst.append(lft_bound)
                    x_to_lst.append(rgt_bound)
                    y_from_lst.append(bot_bound)
                    y_to_lst.append(top_bound)

    return x_from_lst, x_to_lst, y_from_lst, y_to_lst

def find_bound_3d(diffdata3d):
    n_frm = numpy.size(diffdata3d[:, 0:1, 0:1])
    n_row = numpy.size(diffdata3d[0:1, :, 0:1])
    n_col = numpy.size(diffdata3d[0:1, 0:1, :])

    tmp_3d_mask = numpy.zeros(n_row * n_col * n_frm , dtype = int).reshape(n_frm, n_row, n_col)
    tmp_3d_mask[:, :, :] = diffdata3d[:, :, :]
    x_from_lst = []
    x_to_lst = []
    y_from_lst = []
    y_to_lst = []
    z_from_lst = []
    z_to_lst = []
    for frm in range(0, n_frm, 1):
        for row in range(0, n_row, 1):
            for col in range(0, n_col, 1):
                if tmp_3d_mask[frm, row, col] == 1:
                    bck_bound = frm - 1
                    frn_bound = frm + 1
                    lft_bound = col - 1
                    rgt_bound = col + 1
                    bot_bound = row - 1
                    top_bound = row + 1

                    bck_old = 0
                    frn_old = 0
                    bot_old = 0
                    top_old = 0
                    lft_old = 0
                    rgt_old = 0
                    in_img = "True"
                    while in_img == "True":

                        # left wall
                        if lft_bound >= 0 and top_bound < n_row and bot_bound >= 0 and frn_bound < n_frm and bck_bound >= 0:
                            stay_in = "True"
                            for scan_row in range(bot_bound, top_bound + 1, 1):
                                for scan_frm in range(bck_bound, frn_bound + 1, 1):
                                    if tmp_3d_mask[scan_frm, scan_row, lft_bound] == 1:
                                        lft_bound -= 1
                                        stay_in = "False"
                                        break
                                if stay_in == "False":
                                    break
                        else:
                            in_img = "false"
                            break

                        # right wall
                        if rgt_bound < n_col and top_bound < n_row and bot_bound >= 0 and frn_bound < n_frm and bck_bound >= 0:
                            stay_in = "True"
                            for scan_row in range(bot_bound, top_bound + 1, 1):
                                for scan_frm in range(bck_bound, frn_bound + 1, 1):
                                    if tmp_3d_mask[scan_frm, scan_row, rgt_bound] == 1:
                                        rgt_bound += 1
                                        stay_in = "False"
                                        break
                                if stay_in == "False":
                                    break
                        else:
                            in_img = "false"
                            break

                        # bottom wall
                        if bot_bound >= 0 and lft_bound >= 0 and rgt_bound < n_col and frn_bound < n_frm and bck_bound >= 0:
                            stay_in = "True"
                            for scan_col in range(lft_bound, rgt_bound + 1 , 1):
                                for scan_frm in range(bck_bound, frn_bound + 1, 1):
                                    if tmp_3d_mask[scan_frm, bot_bound, scan_col] == 1:
                                        bot_bound -= 1
                                        stay_in = "False"
                                        break
                                if stay_in == "False":
                                    break
                        else:
                            in_img = "false"
                            break
                        # top wall
                        if top_bound < n_row and lft_bound >= 0 and rgt_bound < n_col and frn_bound < n_frm and bck_bound >= 0:
                            stay_in = "True"
                            for scan_col in range(lft_bound, rgt_bound + 1 , 1):
                                for scan_frm in range(bck_bound, frn_bound + 1, 1):
                                    if tmp_3d_mask[scan_frm, top_bound, scan_col] == 1:
                                        top_bound += 1
                                        stay_in = "False"
                                        break
                                if stay_in == "False":
                                    break
                        else:
                            in_img = "false"
                            break

                        # front wall
                        if frn_bound < n_frm and lft_bound >= 0 and rgt_bound < n_col and top_bound < n_row and bot_bound >= 0:
                            stay_in = "True"
                            for scan_col in range(lft_bound, rgt_bound + 1 , 1):
                                for scan_row in range(bot_bound, top_bound + 1, 1):
                                    if tmp_3d_mask[frn_bound, scan_row, scan_col] == 1:
                                        frn_bound += 1
                                        stay_in = "False"
                                        break
                                if stay_in == "False":
                                    break
                        else:
                            in_img = "false"
                            break

                        # back wall
                        if bck_bound >= 0 and lft_bound >= 0 and rgt_bound < n_col and top_bound < n_row and bot_bound >= 0:
                            stay_in = "True"
                            for scan_col in range(lft_bound, rgt_bound + 1 , 1):
                                for scan_row in range(bot_bound, top_bound + 1, 1):
                                    if tmp_3d_mask[bck_bound, scan_row, scan_col] == 1:
                                        bck_bound -= 1
                                        stay_in = "False"
                                        break
                                if stay_in == "False":
                                    break
                        else:
                            in_img = "false"
                            break

                        if bot_bound == bot_old and top_bound == top_old and lft_bound == lft_old and rgt_bound == rgt_old and bck_bound == bck_old and frn_bound == frn_old:
                            break
                        bot_old = bot_bound
                        top_old = top_bound
                        lft_old = lft_bound
                        rgt_old = rgt_bound
                        bck_old = bck_bound
                        frn_old = frn_bound
                    #lst_coord.append([bot_bound, lft_bound, top_bound, rgt_bound])

                    if  in_img == "True":
                        tmp_3d_mask[bck_bound:frn_bound, bot_bound:top_bound, lft_bound:rgt_bound] = 0
                        x_from_lst.append(lft_bound)
                        x_to_lst.append(rgt_bound)
                        y_from_lst.append(bot_bound)
                        y_to_lst.append(top_bound)
                        z_from_lst.append(bck_bound)
                        z_to_lst.append(frn_bound)

    return x_from_lst, x_to_lst, y_from_lst, y_to_lst, z_from_lst, z_to_lst
