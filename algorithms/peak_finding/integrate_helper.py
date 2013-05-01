from scitbx.array_family import flex
from dials.algorithms.peak_finding import ref_2d
import numpy

n_times = 5
dat2d = numpy.zeros((10, 10), dtype = int)
print dat2d
ref_2d()
# dat2d_ref[:, :] = ref_2d(flex.int(dat2d), n_times).as_numpy_array()
print dat2d
