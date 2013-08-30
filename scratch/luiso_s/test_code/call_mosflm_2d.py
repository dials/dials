from dials.scratch.luiso_s.test_code.mosflm_2D import calc_background_n_make_2d_profile
from dials.scratch.luiso_s.test_code.mosflm_2D import fit_profile_2d
def mosflm_caller(rlist):

    profile, tr_hold = calc_background_n_make_2d_profile(rlist)

    fit_profile_2d(rlist, profile, tr_hold)

