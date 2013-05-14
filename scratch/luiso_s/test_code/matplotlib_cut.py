from matplotlib import pylab, cm
pylab.subplot2grid((2, 3), (0, 0), colspan=3)
pylab.plot(line)
pylab.subplot2grid((2, 3), (1, 0))
pylab.imshow(grid.as_numpy_array(), interpolation='none')
pylab.subplot2grid((2, 3), (1, 1))
pylab.imshow(mask_normal.as_numpy_array(), interpolation='none')
pylab.subplot2grid((2, 3), (1, 2))
pylab.imshow(mask_poisson.as_numpy_array(), interpolation='none')
pylab.show()

