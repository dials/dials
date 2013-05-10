from scitbx.array_family import flex
from dials.algorithms.peak_finding import ref_2d

import numpy

n_times = 5
dat2d = numpy.zeros((100, 100), dtype = numpy.int32)

dat2d[:, :] = 11
print dat2d

dat2d_flex = flex.int(dat2d)
flex_ref2d = ref_2d(dat2d_flex, 10, 20, .33333, 1000, 0.5)
dat2d_ref = flex_ref2d.as_numpy_array()

from matplotlib import pyplot as plt
plt.imshow(dat2d_ref , interpolation = "nearest")
plt.show()

print dat2d_ref
