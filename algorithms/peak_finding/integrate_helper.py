from scitbx.array_family import flex
from dials.algorithms.peak_finding import ref_2d

import numpy

n_times = 5
dat2d = numpy.zeros((100, 100), dtype = numpy.int32)

dat2d[:, :] = 11
print dat2d

dat2d_flex = flex.int(dat2d)

flex_ref2d_g = ref_2d(dat2d_flex, 10, 20, .33333, 1000, 0.9)
dat2d_ref_g = flex_ref2d_g.as_numpy_array()


flex_ref2d_l = ref_2d(dat2d_flex, 10, 20, .33333, 1000, 0.1)
dat2d_ref_l = flex_ref2d_l.as_numpy_array()



from matplotlib import pyplot as plt

plt.imshow(dat2d_ref_g , interpolation = "nearest")
plt.show()

plt.imshow(dat2d_ref_l , interpolation = "nearest")
plt.show()


