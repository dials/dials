from scitbx.array_family import flex
from dials.algorithms.peak_finding import ref_2d

flex_ref2d = ref_2d(20, 20, 1, 2, .33333, 125, 0.5)

dat2d_ref = flex_ref2d.as_numpy_array()
print dat2d_ref
from matplotlib import pyplot as plt

plt.imshow(dat2d_ref , interpolation = "nearest")
plt.show()


