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

for xpos in range(80):
  for ypos in range(86):
    row_str = ypos * nrow
    col_str = xpos * ncol
    dx = 1200 - col_str
    dy = 1300 - row_str
    ref_ang = float( math.atan2(dy, dx) / 3.14159265358)
    i_loc = random.randint(0,999)
    thold = i_loc/20
    ref2d = model_2d(nrow, ncol, 5, 1, ref_ang, i_loc, 0.5)

    '''
    flex_double model_2d(int nrow, int ncol, float a, float b,
                float delta_ang, float imax, float asp)
    '''
    data2d_tmp = ref2d.as_numpy_array()
    data2d[col_str:col_str + ncol, row_str:row_str + nrow] += numpy.float64(data2d_tmp)

    new_r = Reflection()
    new_r.bounding_box = [col_str, col_str + ncol, row_str, row_str + nrow, 0, 1]
    new_r.centroid_position = [col_str + 14.5, row_str + 14.5, 0.5]
    new_r.image_coord_px = [col_str + 14.5, row_str + 14.5]
    #new_r.intensity = 50
    #new_r.intensity_variance = 70

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

from matplotlib import pyplot as plt
print "Plotting data2d"
plt.imshow(data2d, interpolation = "nearest")
plt.show()



from dials.algorithms.background.inclined_background_subtractor \
  import layering_and_background_plane
layering_and_background_plane(rlist)
from dials.algorithms.integration import flex_2d_layering_n_integrating
from dials.algorithms.integration.call_mosflm_2d  import mosflm_caller
flex_2d_layering_n_integrating(rlist)
rlist = mosflm_caller(rlist, xmax, ymax, 3)