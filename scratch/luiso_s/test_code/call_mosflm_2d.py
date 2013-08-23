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
from dials.scratch.luiso_s import model_2d
import numpy
import random
from dials.model.data import Reflection, ReflectionList
from scitbx.array_family import flex
rlist = ReflectionList()

for ivar in range(16):
    nrow = 16
    ncol = 18
    data2d = numpy.zeros((nrow, ncol), dtype = numpy.float64)

    ref_ang = float((ivar / 30) + .65)
    ref2d = model_2d(12, 12, 3, 5, ref_ang, ivar * 5 + 80 , 0.5) # make the 100 smaller than 50 to make it crash
    data2d[3:15, 3:15] += numpy.float64(ref2d.as_numpy_array())
    data2d[:, :] += 20

    mask2d = numpy.zeros((nrow, ncol), dtype = numpy.int32)

    bkg_const = 30
    background2d = numpy.copy(data2d)
    for row in range(nrow):
        for col in range(ncol):
            background2d[row, col] = bkg_const + random.random() * bkg_const / 3 + row + col


    for row in range(nrow):
        for col in range(ncol):
            if(data2d[row, col] < background2d[row, col]):
                data2d[row, col] = background2d[row, col]
                mask2d[row, col] = 0
            else:
                mask2d[row, col] = 1


    data3d = data2d
    data3d.shape = (1,) + data2d.shape
    #print data3d.shape

    mask3d = mask2d
    mask3d.shape = (1,) + mask2d.shape
    #print mask3d.shape

    background3d = background2d
    background3d.shape = (1,) + background2d.shape
    #print background3d.shape


    r = Reflection()
    r.shoebox = flex.double(data3d)
    r.shoebox_mask = flex.int(mask3d)
    r.shoebox_background = flex.double(background3d)
    r.centroid_position = (9.5, 9.5, 0.5)
    #r.intensity = 150.0

    rlist.append(r)





from dials.algorithms.integration import flex_2d_layering_n_integrating
flex_2d_layering_n_integrating(rlist)

old_rlist = numpy.zeros((len(rlist)), dtype = numpy.float64)
old_rlist_var = numpy.zeros((len(rlist)), dtype = numpy.float64)
for n in range(len(rlist)):
    print rlist[n].intensity, rlist[n].intensity_variance
    old_rlist[n] = rlist[n].intensity
    old_rlist_var[n] = rlist[n].intensity_variance

from dials.scratch.luiso_s.test_code.mosflm_2D import calc_background_n_make_2d_profile
profile, tr_hold = calc_background_n_make_2d_profile(rlist)

#from dials.scratch.luiso_s import write_2d
#print "profile ="
#write_2d(profile)

print "tr_hold =", tr_hold


from dials.scratch.luiso_s.test_code.mosflm_2D import fit_profile_2d
fit_profile_2d(rlist, profile, tr_hold)


print "____________________________________________________________________________"
print "sum integration                      vs                      Profile fitting"
print "____________________________________________________________________________"
print "I             Variance              vs           I            Variance     "
for n in range(len(rlist)):
    print old_rlist[n], ", ", old_rlist_var[n], "       vs         ", rlist[n].intensity, ", ", rlist[n].intensity_variance


from matplotlib import pyplot as plt
for r in rlist:
    data2d = r.shoebox.as_numpy_array()
    matrix_img = numpy.zeros((r.shoebox.all()[1] , r.shoebox.all()[2]), dtype = numpy.float64)
    matrix_img[:, :] = data2d[0:1, :, :]
    plt.imshow(matrix_img, interpolation = "nearest")
    plt.show()
    '''
    data2d = r.shoebox_mask.as_numpy_array()
    matrix_img = numpy.zeros((r.shoebox.all()[1] , r.shoebox.all()[2]), dtype = numpy.float64)
    matrix_img[:, :] = data2d[0:1, :, :]
    plt.imshow(matrix_img, interpolation = "nearest")
    plt.show()
    '''
    data2d = r.shoebox_background.as_numpy_array()
    matrix_img = numpy.zeros((r.shoebox.all()[1] , r.shoebox.all()[2]), dtype = numpy.float64)
    matrix_img[:, :] = data2d[0:1, :, :]
    plt.imshow(matrix_img, interpolation = "nearest")
    plt.show()


data2d_prof = profile.as_numpy_array()
plt.imshow(data2d_prof, interpolation = "nearest")
plt.show()
