from __future__ import division
from dials.model.data import Reflection, ReflectionList
from dials.scratch.luiso_s import add_2d, write_2d, subtrac_bkg_2d, fitting_2d
from dials.algorithms.background import curved_background_flex_2d
from scitbx.array_family import flex


def calc_background_n_make_2d_profile(reflections):
    big_nrow = 0
    big_ncol = 0
    counter = 0
    select_rlist = ReflectionList()
    max_i = 0.0
    for ref in reflections:
        if ref.is_valid():
            if ref.intensity > max_i:
                max_i = ref.intensity
    print "________ max_i =", max_i
    thold = 0.5 * max_i
    print "________ thold =", thold
    for ref in reflections:
        if ref.is_valid() and ref.intensity > thold:
            print "I =", ref.intensity
            select_rlist.append(ref)



    for ref in select_rlist:
        shoebox = ref.shoebox
        mask = ref.shoebox_mask
        background = ref.shoebox_background

        data2d = shoebox[0:1, :, :]
        mask2d = mask[0:1, :, :]

        data2d.reshape(flex.grid(shoebox.all()[1:]))
        mask2d.reshape(flex.grid(shoebox.all()[1:]))

        background2d = curved_background_flex_2d(data2d.as_double(), mask2d)
        background2d.reshape(flex.grid(1, background2d.all()[0], background2d.all()[1]))
        background[0:1, :, :] = background2d.as_double()
        local_nrow = shoebox.all()[1]
        local_ncol = shoebox.all()[2]
        if local_nrow > big_nrow:
            big_nrow = local_nrow
        if local_ncol > big_ncol:
            big_ncol = local_ncol
        counter += 1

    print "big (nrow, ncol) =", big_nrow, big_ncol
    big_nrow = big_nrow * 2 + 1
    big_ncol = big_ncol * 2 + 1
    print "big (nrow, ncol) used for profile =", big_nrow, big_ncol
    sumation = flex.double(flex.grid(big_nrow, big_ncol))
    descr = flex.double(flex.grid(1, 3))
    for ref in select_rlist:
        shoebox = ref.shoebox
        mask = ref.shoebox_mask
        background = ref.shoebox_background
        data2d = shoebox[0:1, :, :]
        mask2d = mask[0:1, :, :]
        background2d = background[0:1, :, :]
        data2d.reshape(flex.grid(shoebox.all()[1:]))
        mask2d.reshape(flex.grid(shoebox.all()[1:]))
        background2d.reshape(flex.grid(shoebox.all()[1:]))

        #print "ref.centroid_position =", ref.centroid_position
        descr[0, 0] = ref.centroid_position[0]
        descr[0, 1] = ref.centroid_position[1]
        descr[0, 2] = 1.0 / (ref.intensity * counter)
        #print "background2d ="
        peak2d = subtrac_bkg_2d(data2d, background2d)
        #print "peak 2d ="

        sumation = add_2d(descr, peak2d, sumation)


    return sumation, thold
def fit_profile_2d(reflections, average, thold):
    descr = flex.double(flex.grid(1, 3))
    for ref in reflections:
        if ref.is_valid() and ref.intensity < thold:
            shoebox = ref.shoebox
            mask = ref.shoebox_mask
            background = ref.shoebox_background

            data2d = shoebox[0:1, :, :]
            mask2d = mask[0:1, :, :]
            background2d = background[0:1, :, :]

            data2d.reshape(flex.grid(shoebox.all()[1:]))
            mask2d.reshape(flex.grid(shoebox.all()[1:]))
            background2d.reshape(flex.grid(shoebox.all()[1:]))

            descr[0, 0] = ref.centroid_position[0]
            descr[0, 1] = ref.centroid_position[1]
            descr[0, 2] = 1.0 #/ (ref.intensity * counter)

            I_R = fitting_2d(descr, data2d, background2d, average)
            ref.intensity = I_R[0]
            ref.intensity_variance = I_R[1]
            #return I_R[0], I_R[1]
    return reflections
