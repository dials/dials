from __future__ import division

from scitbx.array_family import flex
from dials.scratch.luiso_s import model_2d
from matplotlib import pyplot as plt
from dials.scratch.luiso_s import add_2d, write_2d
'''
import numpy

data2d = numpy.zeros((40, 60), dtype = numpy.float64)
nrow = 10
ncol = 10


sumation = flex.double(flex.grid(21, 21))
descr = flex.double(flex.grid(1, 3))
descr[0, 0] = .5
descr[0, 1] = 5
descr[0, 2] = 5
print "____________________________________________________________________"


for xpos in range(3):
    for ypos in range(3):
        row_str = ypos * 12
        col_str = xpos * 20
        ref_ang = float((ypos + 1) / 10)
        #flex_int model_2d(int nrow, int ncol, float a, float b,
        #float delta_ang, float imax, float asp)
        ref2d = model_2d(nrow, ncol, 2, 1, ref_ang, 55, 0.5)
        data2d_tmp = ref2d.as_numpy_array()
        data2d[row_str:row_str + nrow, col_str:col_str + ncol] += numpy.float64(data2d_tmp)

        sumation = add_2d(descr, flex.double (numpy.float64 (data2d_tmp)), sumation)
        write_2d(sumation)


print "Plotting data2d"
plt.imshow(data2d, interpolation = "nearest")
plt.show()

print "Plotting reslt"
img_suma = sumation.as_numpy_array()
plt.imshow(img_suma, interpolation = "nearest")
plt.show()
'''
def calc_background_n_make_2d_profile(reflections):
    print "Hi 01"
    big_nrow = 0
    big_ncol = 0
    from dials.algorithms.background import curved_background_flex_2d
    from scitbx.array_family import flex
    for ref in reflections:
        if ref.is_valid():
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

    print "big (nrow, ncol) =", big_nrow, big_ncol
    big_nrow = big_nrow * 2 + 1
    big_ncol = big_ncol * 2 + 1
    sumation = flex.double(flex.grid(big_nrow, big_ncol))
    descr = flex.double(flex.grid(1, 3))
    for ref in reflections:
        if ref.is_valid():
            shoebox = ref.shoebox
            mask = ref.shoebox_mask
            background = ref.shoebox_background

            data2d = shoebox[0:1, :, :]
            mask2d = mask[0:1, :, :]
            background2d = background[0:1, :, :]
            data2d.reshape(flex.grid(shoebox.all()[1:]))
            mask2d.reshape(flex.grid(shoebox.all()[1:]))
            background2d.reshape(flex.grid(shoebox.all()[1:]))

            #print "ref.centroid =", ref.centroid
            print "ref.centroid_position =", ref.centroid_position
            descr[0, 0] = ref.centroid_position[0]
            descr[0, 1] = ref.centroid_position[1]
            descr[0, 2] = 1.0
            sumation = add_2d(descr, data2d, sumation)

    print "Plotting reslt"
    img_suma = sumation.as_numpy_array()
    plt.imshow(img_suma, interpolation = "nearest")
    plt.show()
    print "hi 02"
