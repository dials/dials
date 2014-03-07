from __future__ import division
from dials.array_family import flex
from dials.scratch.luiso_s import model_2d
from matplotlib import pylab

import numpy
import math
import random

xmax = 2400
ymax = 2600

data2d = numpy.zeros((xmax, ymax), dtype = numpy.float64)
nrow = 30
ncol = 30

pi = 3.14159265358

n_x = 80
n_y = 86
num_ref = n_y * n_x
ref_table = flex.reflection_table()

t_shoebox = flex.shoebox(num_ref)
t_intensity = flex.double(num_ref)
t_intensity_var = flex.double(num_ref)
t_bbox = flex.int6(num_ref)
t_xyzobs = flex.vec3_double(num_ref)
t_xyzcal = flex.vec3_double(num_ref)


t_row = 0
for ypos in range(n_y):
  for xpos in range(n_x):
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
      t_bbox[t_row] = [col_str, col_str + ncol, row_str, row_str + nrow, 0, 1]
      t_xyzobs[t_row] = [col_str + 14.5, row_str + 14.5, 0.5]


      t_xyzcal[t_row] = [col_str + 14.5, row_str + 14.5, 0.5]
      np_shoebox = numpy.copy(data2d[col_str:col_str + ncol, row_str:row_str + nrow])
      fl_shoebox = flex.double(np_shoebox)
      fl_shoebox.reshape(flex.grid(1, ncol, nrow))
      fl_shoebox_bkg=fl_shoebox[:,:,:]
      t_shoebox[t_row].data = fl_shoebox
      t_shoebox[t_row].background = fl_shoebox_bkg
      lc_mask = flex.int(flex.grid(1, ncol, nrow), 3)
      for x_loc in range(ncol):
        for y_loc in range(nrow):
          if ref2d[y_loc,x_loc]>thold:
            lc_mask[0, y_loc, x_loc] = 5
      t_shoebox[t_row].mask = lc_mask
      t_row += 1

print "t_row =", t_row
ref_table['shoebox'] = t_shoebox
ref_table['intensity.raw.value'] = t_intensity
ref_table['intensity.raw.variance'] = t_intensity_var
ref_table['bbox'] = t_bbox
ref_table['xyzobs.px.value'] = t_xyzobs
ref_table['xyzcal.px'] = t_xyzcal


ref2d = model_2d(xmax, ymax, 380, 740, 0.25, 955, 0.5)
data2d_tmp = ref2d.as_numpy_array()
data2d[:, :] += numpy.float64(data2d_tmp)

from matplotlib import pyplot as plt
print "Plotting data2d"
plt.imshow(data2d, interpolation = "nearest")
plt.show()


from dials.algorithms.background.inclined_background_subtractor \
  import layering_and_background_plane
layering_and_background_plane(ref_table)
from dials.algorithms.integration import flex_2d_layering_n_integrating
flex_2d_layering_n_integrating(ref_table)


print "_____________________________________________________ here"

#t_intensity = ref_table['intensity.raw.value']

old_i_table = t_intensity[:]

for tmp_i in (t_intensity):
  print "tmp_i = ", tmp_i
#tmp='''
print "adding noise ...."
t_row = 0
for count in range(num_ref):
    for x_loc in range(ncol):
      for y_loc in range(nrow):
        roll_the_dice = random.randint(1,50)
        if roll_the_dice <=5:
          t_shoebox[t_row].data[0, y_loc, x_loc] = -1
          t_shoebox[t_row].mask[0, y_loc, x_loc] = 0
        else:
          t_shoebox[t_row].data[0, y_loc, x_loc] += random.randint(0,10)

    t_row += 1
print "adding noise .... done"
#'''

layering_and_background_plane(ref_table)
flex_2d_layering_n_integrating(ref_table)

paint_compare = []
for i in range(len(t_intensity)):
  paint_compare.append([ t_intensity[i], old_i_table[i]])
paint_compare_sort = sorted(paint_compare)

data1d = numpy.zeros(num_ref, dtype = numpy.float64)
new_data1d = numpy.zeros(num_ref, dtype = numpy.float64)
for i in range(num_ref):
  data1d[i] = paint_compare_sort[i][0]
  new_data1d[i] = paint_compare_sort[i][1]

from matplotlib import pylab
pylab.plot(data1d)
pylab.plot(new_data1d)
pylab.show()

temp = '''
old_r_list = rlist[:]

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
  paint_compare.append([rlist[i].intensity, rlist[i].intensity_variance, \
              old_r_list[i].intensity, old_r_list[i].intensity_variance])

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
###########################################################################
paint_compare = []
paint_compare_sort = []
for i in range(len(rlist)):
  paint_compare.append([ rlist[i].intensity, old_r_list[i].intensity])

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
'''