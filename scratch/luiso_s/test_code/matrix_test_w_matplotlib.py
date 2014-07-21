from __future__ import division
import numpy

ysize = 10
xsize = 16
mat2d = numpy.arange( ysize * xsize, dtype = 'uintc' ).reshape( ysize, xsize )

x = 5
mat2d[0:5, 2:6] = 50

coords = [3.5,1.5]

import matplotlib.pyplot as plt
plt.imshow(numpy.transpose(mat2d), interpolation = "nearest")
ini_x = ini_y = - 0.5
end_y = xsize - 0.5
end_x = ysize - 0.5

coords[0] = coords[0] - 0.5
coords[1] = coords[1] - 0.5

x_lin_from = (ini_x + coords[0]) / 2.0
y_lin_from = (ini_y + coords[1]) / 2.0

x_lin_to = (coords[0] + end_x) / 2.0
y_lin_to = (coords[1] + end_y) / 2.0

plt.vlines(coords[0], y_lin_from, y_lin_to)
plt.hlines(coords[1], x_lin_from, x_lin_to)

#plt.vlines(coords[0] - 0.5, ini_x, end_x)
#plt.hlines(coords[1] - 0.5, ini_y, end_y)

plt.show()
plt.close()