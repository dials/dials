from dials.scratch.luiso_s.test_code.mosflm_2D import calc_background_n_make_2d_profile
from dials.scratch.luiso_s.test_code.mosflm_2D import fit_profile_2d
from dials.model.data import Reflection, ReflectionList
def mosflm_caller(rlist, xmax, ymax, n_div):

    ncol = n_div
    nrow = n_div
    arr_rlist = []
    for col in range(ncol):
        b = []
        for row in range(nrow):
            b.append(ReflectionList())
        arr_rlist.append(b)

    for r in rlist:
        x, y = r.image_coord_px
        col = int(float(x) / float(xmax) * n_div)
        row = int(float(y) / float(ymax) * n_div)
        arr_rlist[row][col].append(r)

    for col in range(ncol):
        for row in range(nrow):
            profile, tr_hold = calc_background_n_make_2d_profile(arr_rlist[row][col])
            arr_rlist[row][col] = fit_profile_2d(arr_rlist[row][col], profile, tr_hold)

    new_rlist = ReflectionList()
    for col in range(ncol):
        for row in range(nrow):
            for r in arr_rlist[row][col]:
                new_rlist.append(r)

    for r in new_rlist:
        print r.intensity, r.intensity_variance

    return new_rlist
