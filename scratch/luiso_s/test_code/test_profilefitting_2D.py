from __future__ import division
from dials.array_family import flex
from dials.scratch.luiso_s import model_2d
from matplotlib import pylab

import numpy
import math
import random

xmax = 240
ymax = 250

data2d = numpy.zeros((xmax, ymax), dtype = numpy.float64)

n_x = 23
n_y = 24

nrow = int(float(ymax) / float(n_y))

ncol = int(float(xmax) / float(n_x))

centr_col = float(ncol)/2.0
centr_row = float(ncol)/2.0

print "ncol =", ncol
print "nrow =", nrow

if ncol > nrow:
  ncol = nrow
elif nrow > ncol:
  nrow = ncol

print "ncol =", ncol
print "nrow =", nrow

pi = 3.14159265358

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
      dx = col_str - xmax / 2
      dy = row_str - ymax / 2

      ref_ang = float( math.atan2(dx, dy) / pi)
      i_loc = random.randint(0,999)
      thold = i_loc / 8
      ref2d = model_2d(nrow, ncol, 5, 1, ref_ang, i_loc, 0.5)

      data2d_tmp = ref2d.as_numpy_array()
      data2d[col_str:col_str + ncol, row_str:row_str + nrow] += \
       numpy.float64(data2d_tmp)
      t_bbox[t_row] = [col_str, col_str + ncol, row_str, row_str + nrow, 0, 1]
      t_xyzobs[t_row] = [col_str + centr_col, row_str + centr_row, 0.5]

      t_xyzcal[t_row] = [col_str + centr_col, row_str + centr_row, 0.5]
      np_shoebox = numpy.copy( \
       data2d[col_str:col_str + ncol, row_str:row_str + nrow])

      fl_shoebox = flex.double(np_shoebox)
      fl_shoebox.reshape(flex.grid(1, ncol, nrow))
      fl_shoebox_bkg=fl_shoebox[:,:,:]
      t_shoebox[t_row].data = fl_shoebox
      t_shoebox[t_row].background = fl_shoebox_bkg
      lc_mask = flex.int(flex.grid(1, ncol, nrow), 3)
      for x_loc in range(ncol):
        for y_loc in range(nrow):
          if ref2d[y_loc,x_loc] > thold:
            lc_mask[0, y_loc, x_loc] = 5
      t_shoebox[t_row].mask = lc_mask
      t_row += 1

print "t_row =", t_row
ref_table['shoebox'] = t_shoebox
ref_table['intensity.raw.value'] = t_intensity
ref_table['intensity.raw.variance'] = t_intensity_var
ref_table['bbox'] = t_bbox
#ref_table['xyzobs.px.value'] = t_xyzobs
ref_table['xyzcal.px'] = t_xyzcal

ref2d = model_2d(xmax, ymax, 380, 740, 0.25, 955, 0.5)
data2d_tmp = ref2d.as_numpy_array()
data2d[:, :] += numpy.float64(data2d_tmp)

#tmp = '''
from matplotlib import pyplot as plt
print "Plotting data2d"
plt.imshow(data2d, interpolation = "nearest")
plt.show()
#'''

from dials.algorithms.background.inclined_background_subtractor \
  import layering_and_background_plane
layering_and_background_plane(ref_table)
from dials.algorithms.integration import flex_2d_layering_n_integrating
flex_2d_layering_n_integrating(ref_table)


t_intensity = ref_table['intensity.raw.value']
old_i_table = t_intensity[:]

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

from dials.algorithms.integration.call_mosflm_2d  import mosflm_caller
pf_ref_table = mosflm_caller(ref_table, xmax, ymax, 3)

t_intensity = pf_ref_table['intensity.raw.value']

paint_compare = []
for i in range(len(t_intensity)):
  paint_compare.append([ old_i_table[i], t_intensity[i]])
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
