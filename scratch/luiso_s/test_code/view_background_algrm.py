from __future__ import division
from scitbx.array_family import flex
from dials.scratch.luiso_s import model_2d

nrow = 35
ncol = 35
ref_ang = 5
local_i = 170

data2d = model_2d(nrow, ncol, 6, 3, ref_ang, local_i, 0.5)

ext_siz = 2
mask2d = flex.int(flex.grid(nrow, ncol))
for x_loc in range(ncol):
  for y_loc in range(nrow):
    if( data2d[y_loc, x_loc] > 10):
      mask2d[y_loc, x_loc] = 5
    else:
      mask2d[y_loc, x_loc] = 3

ext_siz = 1
tmp_mask = mask2d[:,:]
mask2d = mask2d.as_numpy_array()
for x_loc in range(ncol):
  for y_loc in range(nrow):
    if( tmp_mask[y_loc, x_loc] == 5):
      mask2d[y_loc - ext_siz:y_loc + ext_siz \
      , x_loc - ext_siz : x_loc + ext_siz] = 5

mask2d = flex.int(mask2d)
print "adding noise ...."
import random
for x_loc in range(ncol):
  for y_loc in range(nrow):
    roll_the_dice = random.randint(1,50)

    noise_scale = .3

    dx = ((x_loc + y_loc) - ((ncol + nrow)/2)) * 0.3
    parabl = dx * dx
    i_func = parabl * 0.03
    data2d[y_loc, x_loc] += parabl
    data2d[y_loc, x_loc] += roll_the_dice * noise_scale

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
