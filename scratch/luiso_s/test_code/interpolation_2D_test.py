from __future__ import division
from scitbx.array_family import flex
from dials.scratch.luiso_s import model_2d
from dials.model.data import Reflection, ReflectionList
from matplotlib import pylab

import numpy
import math
import random

xmax = 2400
ymax = 2600

data2d = numpy.zeros((xmax, ymax), dtype = numpy.float64)
nrow = 30
ncol = 30

rlist = ReflectionList()
pi = 3.14159265358

for ypos in range(86):
  for xpos in range(80):
    #if xpos/12.0 == int(xpos/12) and ypos/12.0 == int(ypos/12):
    #if ypos/6.0 == int(ypos/6):
    if xpos/3.0 == int(xpos/3) and ypos/3.0 == int(ypos/3):
      row_str = ypos * nrow
      col_str = xpos * ncol
      dx = col_str - 1200
      dy = row_str - 1300

      if dx < 0 and dy < 0:
        ref_ang = 0.25
      elif dx > 0 and dy < 0:
        ref_ang = 0.75
      elif dx < 0 and dy > 0:
        ref_ang = 0.75
      else:
        ref_ang = 0.20
      ref_ang = float( math.atan2(dx, dy) / pi)
      i_loc = random.randint(0,999)
      thold = i_loc/20
      ref2d = model_2d(nrow, ncol, 5, 1, ref_ang, i_loc, 0.5)

      data2d_tmp = ref2d.as_numpy_array()
      data2d[col_str:col_str + ncol, row_str:row_str + nrow] += numpy.float64(data2d_tmp)

      new_r = Reflection()
      new_r.bounding_box = [col_str, col_str + ncol, row_str, row_str + nrow, 0, 1]
      new_r.centroid_position = [col_str + 14.5, row_str + 14.5, 0.5]
      new_r.image_coord_px = [col_str + 14.5, row_str + 14.5]

      np_shoebox = numpy.copy(data2d[col_str:col_str + ncol, row_str:row_str + nrow])
      fl_shoebox = flex.double(np_shoebox)
      fl_shoebox.reshape(flex.grid(1, ncol, nrow))
      fl_shoebox_bkg=fl_shoebox[:,:,:]

      new_r.shoebox = fl_shoebox
      new_r.shoebox_background = fl_shoebox_bkg
      mask = flex.int(flex.grid(1, ncol, nrow), 3)
      for x_loc in range(ncol):
        for y_loc in range(nrow):
          if ref2d[y_loc,x_loc]>thold:
            mask[0, y_loc, x_loc] = 5
      new_r.shoebox_mask = mask

      rlist.append(new_r)

ref2d = model_2d(xmax, ymax, 380, 740, 0.25, 955, 0.5)
data2d_tmp = ref2d.as_numpy_array()
data2d[:, :] += numpy.float64(data2d_tmp)

print "created ", len(rlist), "reflections"

from matplotlib import pyplot as plt
print "Plotting data2d"
plt.imshow(data2d, interpolation = "nearest")
plt.show()


from dials.algorithms.background.inclined_background_subtractor \
  import layering_and_background_plane
layering_and_background_plane(rlist)
from dials.algorithms.integration import flex_2d_layering_n_integrating
flex_2d_layering_n_integrating(rlist)


old_r_list = rlist[:]
tmp='''
print "adding noise ...."
for r in rlist:
    for x_loc in range(ncol):
      for y_loc in range(nrow):
        roll_the_dice = random.randint(1,50)
        if roll_the_dice <=5:
          r.shoebox[0, y_loc, x_loc] = -1
          r.shoebox_mask[0, y_loc, x_loc] = 0
        else:
          r.shoebox[0, y_loc, x_loc] += random.randint(0,10)
print "adding noise .... done"
#'''
from dials.algorithms.background.inclined_background_subtractor \
  import layering_and_background_plane
layering_and_background_plane(rlist)
from dials.algorithms.integration import flex_2d_layering_n_integrating
flex_2d_layering_n_integrating(rlist)
from dials.algorithms.integration.call_mosflm_2d  import mosflm_caller
rlist = mosflm_caller(rlist, xmax, ymax, 3)
paint_compare = []
for i in range(len(rlist)):
  #paint_compare.append([ rlist[i].intensity, old_r_list[i].intensity])
  paint_compare.append([old_r_list[i].intensity, rlist[i].intensity])

paint_compare_sort = sorted(paint_compare)
from matplotlib import pylab
import numpy
data1d = numpy.zeros(len(rlist), dtype = numpy.float64)
data1d_var = numpy.zeros(len(rlist), dtype = numpy.float64)
new_data1d = numpy.zeros(len(rlist), dtype = numpy.float64)
for i in range(len(rlist)):
  data1d[i] = paint_compare_sort[i][0]
  try:
    new_data1d[i] = paint_compare_sort[i][1] / paint_compare_sort[i][0]
  except:
    print ">>>>>>>>>>", paint_compare_sort[i][1] , paint_compare_sort[i][0]
  data1d_var[i] = paint_compare_sort[i][1]

pylab.plot(data1d)
pylab.plot(new_data1d)
pylab.plot(data1d_var)
pylab.show()
