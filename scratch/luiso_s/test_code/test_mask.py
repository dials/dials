from __future__ import division
from scitbx.array_family import flex
from dials.scratch.luiso_s import model_2d


from matplotlib import pylab

import numpy

np_data2d = numpy.zeros((550, 950), dtype = numpy.float64)
nrow = 60
ncol = 60

for xpos in range(20):
  for ypos in range(20):
    row_str = ypos * 25
    col_str = xpos * 45
    ref_ang = float(ypos / 10)
    local_i = (39 - xpos - ypos) * 5
    print "xpos =", xpos
    print "ypos =", ypos
    print "local_i =", local_i
    print "____________________________"
    ref2d = model_2d(nrow, ncol, 3, 2, ref_ang, local_i, 0.5)
    data2d_tmp = ref2d.as_numpy_array()
    np_data2d[row_str:row_str + 60, col_str:col_str + 60] += numpy.float64(data2d_tmp)

ref2d = model_2d(550, 950, 280, 140, 1.0, 955, 0.5)
data2d_tmp = ref2d.as_numpy_array()
np_data2d[:, :] += numpy.float64(data2d_tmp)


data2d = flex.double(np_data2d)
#tmp='''
print "adding noise ...."
import random
for x_loc in range(950):
  for y_loc in range(550):
    roll_the_dice = random.randint(1,50)
    tmp='''
    if roll_the_dice <=2:
      r.shoebox[0, y_loc, x_loc] = -1
      r.shoebox_mask[0, y_loc, x_loc] = 0
    else:
    #'''
    data2d[y_loc, x_loc] += (y_loc + x_loc) * 0.3
    data2d[y_loc, x_loc] += roll_the_dice * 0.0005 * data2d[y_loc, x_loc]
print "adding noise .... done"

from matplotlib import pyplot as plt
print "Plotting data2d"
plt.imshow(data2d.as_numpy_array(), interpolation = "nearest")
plt.show()

#'''
#    code that will become production code:
#    from data2d flex array that contains an image
#    it should return a flex array with the mask

from dials.algorithms.peak_finding import smooth_2d
from dials.algorithms.peak_finding import find_mask_2d
n_times = 5
data2dsmoth = smooth_2d(data2d, n_times)

mask2d = find_mask_2d(data2d, data2dsmoth, n_times)


# end code that will become production code
ploting_smooth = '''
print "Plotting data2dsmoth"
np_data2dsmoth = data2dsmoth.as_numpy_array()
plt.imshow(np_data2dsmoth, interpolation = "nearest")#, cmap = pylab.gray())
plt.show()
'''
print "Plotting data2d mask"
np_data2dmask = mask2d.as_numpy_array()
plt.imshow(np_data2dmask, interpolation = "nearest", cmap = pylab.gray())#, cmap = pylab.gray())
plt.show()