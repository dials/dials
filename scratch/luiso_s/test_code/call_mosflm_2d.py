from __future__ import division
#
#  caller of  mosflm_2d testing code
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
import numpy
'''
from scitbx.array_family import flex

from dials.scratch.luiso_s import add_2d, write_2d

#data2d = numpy.zeros((40, 60), dtype = numpy.float64)
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
        #data2d[row_str:row_str + nrow, col_str:col_str + ncol] += numpy.float64(data2d_tmp)

        sumation = add_2d(descr, flex.double (numpy.float64 (data2d_tmp)), sumation)

'''

#  numpy_integr_2d testing code
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from dials.model.data import Reflection, ReflectionList
from scitbx.array_family import flex
rlist = ReflectionList()

for ivar in range(3):
    print ivar
    nrow = 5 + ivar
    ncol = 7 - ivar
    print nrow, ncol


    data2d = numpy.zeros((nrow, ncol), dtype = numpy.float64)
    data2d[0, 1] = data2d[1, 0] = 1
    data2d[0, 2] = data2d[2, 0] = 3
    data2d[0, 3] = data2d[3, 0] = 2
    data2d[0, 4] = data2d[4, 0] = 4
    data2d[4, 1] = data2d[1, 4] = 1
    data2d[4, 2] = data2d[2, 4] = 3
    data2d[4, 3] = data2d[3, 4] = 2
    data2d[0, 0] = data2d[4, 4] = 5

    data2d[1:4, 1:4] = 10
    data2d[2:3, 2:3] = 50
    print data2d

    mask2d = numpy.zeros((nrow, ncol), dtype = numpy.int32)
    mask2d[1:4, 1:4] = 1
    print mask2d

    background2d = numpy.copy(data2d)
    background2d[:, :] = 0.0

    data3d = data2d
    data3d.shape = (1,) + data2d.shape
    print data3d.shape

    mask3d = mask2d
    mask3d.shape = (1,) + mask2d.shape
    print mask3d.shape

    background3d = background2d
    background3d.shape = (1,) + background2d.shape
    print background3d.shape


    r = Reflection()
    r.shoebox = flex.double(data2d)
    r.shoebox_mask = flex.int(mask2d)
    r.shoebox_background = flex.double(background2d)


    #r.centroid_position = (0, 2, 2)

    rlist.append(r)

from dials.algorithms.centroid.toy_centroid_Lui import toy_centroid_lui
centroid = toy_centroid_lui(rlist)
rlist = centroid.get_reflections()

from dials.scratch.luiso_s.test_code.mosflm_2D import calc_background_n_make_2d_profile
calc_background_n_make_2d_profile(rlist)



for r in rlist:
    print r
    matrix_img = r.shoebox.as_numpy_array()
    print
    print "shoebox"
    print matrix_img

    matrix_bkg = r.shoebox_background.as_numpy_array()
    print
    print "background"
    print matrix_bkg

    matrix_mask = r.shoebox_mask.as_numpy_array()
    print
    print "mask"
    print matrix_mask
