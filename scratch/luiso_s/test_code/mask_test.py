from __future__ import division
from dials.algorithms.shoebox import build_mask
from scitbx.array_family import flex

img_slice_2d = flex.double(flex.grid(23, 17), 0)




#img_slice = img_3d[batch[i] - 1:batch[i], y_from:y_to, x_from:x_to]
#new_r.shoebox = img_slice.as_double()
#img_slice_2d = img_slice
#img_slice_2d.reshape(flex.grid(img_slice.all()[1], img_slice.all()[2]))
#tmp_2d_mask = build_mask(nx[i], ny[i], nrx[i], nry[i], nc[i], img_slice_2d.as_double())
x = build_mask(23   , 17   , 3     , 2     , 8    , img_slice_2d)
#import numpy
print "mask =", x.as_numpy_array()

'''
show_mask = x.as_numpy_array()
for col in range(23):
    for row in range(17):
        if show_mask[row, col] == 5:
            show_mask[row, col] = 8
        elif show_mask[row, col] == 3:
            show_mask[row, col] = 1
print show_mask
'''
