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

xmax = 500
ymax = 400
rlist = ReflectionList()
num_of_ref = 144
for ivar in range(num_of_ref):
  nrow = 16
  ncol = 18
  data2d = numpy.zeros((nrow, ncol), dtype = numpy.float64)

  ref_ang = float((ivar / 30) + .65)
  ref2d = model_2d(12, 12, 3, 5, ref_ang, ivar * 1.5 + 80 , 0.5) # make the 100 smaller than 50 to make it crash
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
  r.image_coord_px = random.random() * xmax, random.random() * ymax
  rlist.append(r)


from dials.algorithms.integration import flex_2d_layering_n_integrating
flex_2d_layering_n_integrating(rlist)

old_rlist = rlist

from dials.scratch.luiso_s.test_code.call_mosflm_2d import mosflm_caller
rlist = mosflm_caller(rlist, xmax, ymax, 4)

print '_____________________________________________________________________________________________'
print "== new intensity , new intensity variance  //// old  intensity ,  old intensity variance  ======"
for pos in range(len(rlist)):
  print '   ', rlist[pos].intensity, ', ' , rlist[pos].intensity_variance, '              ', old_rlist[pos].intensity, ', ' , old_rlist[pos].intensity_variance

print '_____________________________________________________________________________________________'
print "== new coordinate px                ////               old coordinate px   ======"
for pos in range(len(rlist)):
  print '   ', rlist[pos].image_coord_px, '              ', old_rlist[pos].image_coord_px
