from scitbx.array_family import flex
from dials.algorithms.peak_finding import ref_2d, measure_2d

flex_ref2d = ref_2d(20, 20, 1, 2, .33333, 125, 0.5)
a = 0.0
b = 0.0
c = measure_2d(flex_ref2d, a, b)
print "a, b, c =", a, b, c


dat2d_ref = flex_ref2d.as_numpy_array()
print dat2d_ref
from matplotlib import pyplot as plt

plt.imshow(dat2d_ref , interpolation = "nearest")
plt.show()


