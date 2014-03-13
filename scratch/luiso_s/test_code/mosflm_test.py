from __future__ import division

from scitbx.array_family import flex
from dials.scratch.luiso_s import model_2d, write_2d
from matplotlib import pyplot as plt

# from dials.scratch.luiso_s import add_2d

from dials.algorithms.integration import add_2d

import numpy

data2d = numpy.zeros((40, 60), dtype = numpy.float64)
nrow = 10
ncol = 10

sumation = flex.double(flex.grid(21, 21))
descr = flex.double(flex.grid(1, 3))
descr[0, 0] = .5
descr[0, 1] = 5
descr[0, 2] = 5
print "____________________________________________________________________"

for xpos in range(3):
  for ypos in range(3):
    row_str = ypos * 12
    col_str = xpos * 20
    ref_ang = float((ypos + 1) / 10)
    # flex_int model_2d(int nrow, int ncol, float a, float b,
    # float delta_ang, float imax, float asp)
    ref2d = model_2d(nrow, ncol, 2, 1, ref_ang, 55, 0.5)
    data2d_tmp = ref2d.as_numpy_array()
    data2d[row_str:row_str + nrow, col_str:col_str + ncol] += numpy.float64(data2d_tmp)

    sumation = add_2d(descr, ref2d, sumation)
    write_2d(sumation)

print "Plotting data2d"
plt.imshow(data2d, interpolation = "nearest")
plt.show()

print "Plotting reslt"
img_suma = sumation.as_numpy_array()
plt.imshow(img_suma, interpolation = "nearest")
plt.show()
