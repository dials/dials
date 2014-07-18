from __future__ import division
import numpy

ysize = xsize = 16
mat2d = numpy.arange( ysize * xsize, dtype = 'uintc' ).reshape( ysize, xsize )

x = 5
mat2d[0:5, 2:6] = 50

import matplotlib.pyplot as plt
plt.imshow(numpy.transpose(mat2d), interpolation = "nearest")
plt.vlines(1, 1, 8)
plt.show()
plt.close()