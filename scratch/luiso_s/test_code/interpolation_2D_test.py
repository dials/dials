from __future__ import division
from scitbx.array_family import flex
from dials.scratch.luiso_s import model_2d
from dials.model.data import Reflection, ReflectionList
from matplotlib import pylab

import numpy
import math
import random
data2d = numpy.zeros((2400, 2600), dtype = numpy.float64)
nrow = 30
ncol = 30

rlist = ReflectionList()

for xpos in range(86):
  for ypos in range(80):
    row_str = ypos * 30
    col_str = xpos * 30
    dx = 1300 - col_str
    dy = 1200 - row_str
    ref_ang = float( math.atan2(dy, dx) / 3.14159265358)
    i_loc = random.randint(0,999)
    thold = i_loc/20
    ref2d = model_2d(nrow, ncol, 5, 1, ref_ang, i_loc, 0.5)

    '''

    flex_double model_2d(int nrow, int ncol, float a, float b,
                float delta_ang, float imax, float asp)
    '''
    data2d_tmp = ref2d.as_numpy_array()
    data2d[row_str:row_str + 30, col_str:col_str + 30] += numpy.float64(data2d_tmp)

    new_r = Reflection()
    new_r.bounding_box = [col_str, col_str + 30, row_str, row_str + 30, 0, 1]
    new_r.centroid_position = [col_str + 14.5, row_str + 14.5, 0.5]
    new_r.intensity = 50
    new_r.intensity_variance = 70



    np_shoebox = numpy.copy(data2d[row_str:row_str + 30, col_str:col_str + 30])
    fl_shoebox = flex.double(np_shoebox)
    fl_shoebox.reshape(flex.grid(1, 30, 30))
    fl_shoebox_bkg=fl_shoebox[:,:,:]

    new_r.shoebox = fl_shoebox
    new_r.shoebox_background = fl_shoebox_bkg
    mask = flex.int(flex.grid(1, 30, 30), 3)
    for x_loc in range(30):
      for y_loc in range(30):
        if ref2d[y_loc,x_loc]>thold:
          mask[0, y_loc, x_loc] = 5
    new_r.shoebox_mask = mask


    from matplotlib import pyplot as plt
    plt.imshow(np_shoebox, interpolation = "nearest")
    plt.show()

    tmp_2d_mask = mask[:, :, :]
    tmp_2d_mask.reshape(flex.grid(30, 30))
    np_2d_mask = tmp_2d_mask.as_numpy_array()
    plt.imshow(np_2d_mask)
    plt.show()
    #img_slice = flex.double(flex.grid(1, y1 - y0, x1 - x0), 0)


    '''
    new_r.shoebox_background = img_slice.as_double()

    img_slice_2d = img_slice
    img_slice_2d.reshape(flex.grid(img_slice.all()[1], img_slice.all()[2]))
    tmp_2d_mask = build_mask(nx[j], ny[j], nrx[j], nry[j], nc[j], img_slice_2d.as_double())

    block_col = int((float(xdet[i]) / float(ncol)) * 5.0 + 1)
    block_row = int((float(ydet[i]) / float(nrow)) * 5.0 + 1)
    box = block_row + (block_col - 1) * 5


    tmp_2d_mask.reshape(flex.grid(1, tmp_2d_mask.all()[0], tmp_2d_mask.all()[1]))
    tmp_2d_mask[0:1, :, :] = tmp_2d_mask
    new_r.shoebox_mask = tmp_2d_mask
    '''
    rlist.append(new_r)

ref2d = model_2d(2400, 2600, 380, 740, 0.25, 955, 0.5)
data2d_tmp = ref2d.as_numpy_array()
data2d[:, :] += numpy.float64(data2d_tmp)

from matplotlib import pyplot as plt
print "Plotting data2d"
plt.imshow(data2d, interpolation = "nearest")
plt.show()





that_should_run = '''
from dials.algorithms.integration import flex_2d_layering_n_integrating
from dials.algorithms.integration.call_mosflm_2d  import mosflm_caller
flex_2d_layering_n_integrating(rlist)
xmax, ymax = sweep.get_detector()[0].get_image_size()
rlist = mosflm_caller(rlist, xmax, ymax, self.nblocks)
'''