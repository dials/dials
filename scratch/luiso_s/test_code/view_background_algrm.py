from __future__ import division
from scitbx.array_family import flex
from dials.scratch.luiso_s import model_2d

from matplotlib import pylab

import numpy

nrow = 40
ncol = 40

ref_ang = 5
local_i = 200

data2d = model_2d(nrow, ncol, 3, 2, ref_ang, local_i, 0.5)

mask2d = flex.int(flex.grid(nrow, ncol))
for x_loc in range(ncol):
  for y_loc in range(nrow):
    if( data2d[y_loc, x_loc] > 10):
      mask2d[y_loc, x_loc] = 5
    else:
      mask2d[y_loc, x_loc] = 3

#tmp='''
print "adding noise ...."
import random
for x_loc in range(ncol):
  for y_loc in range(nrow):
    roll_the_dice = random.randint(1,50)
    data2d[y_loc, x_loc] += (y_loc + x_loc) * 0.3
    noise_scale = .01
    inclined_pl_scale = 30
    data2d[y_loc, x_loc] += inclined_pl_scale * (nrow + ncol)
    data2d[y_loc, x_loc] += roll_the_dice * noise_scale #* data2d[y_loc, x_loc]

    
    
print "adding noise .... done"

from dials.algorithms.background import curved_background_flex_2d
background2d = curved_background_flex_2d(data2d, mask2d)

from matplotlib import pyplot as plt
data2d = data2d.as_numpy_array()
plt.imshow(data2d, interpolation = "nearest")
plt.show()
np_data2dmask = mask2d.as_numpy_array()
plt.imshow(np_data2dmask, interpolation = "nearest")
plt.show()
background2d = background2d.as_numpy_array()
plt.imshow(background2d, interpolation = "nearest")
plt.show()