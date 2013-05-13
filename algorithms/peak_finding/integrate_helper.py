from scitbx.array_family import flex
from dials.algorithms.peak_finding import model_2d, measure_2d_angl

ref2d = model_2d(20, 20, 1, 2, .33333, 125, 0.5)

ang = measure_2d_angl(ref2d, 10, 10)

print "ang =", ang

dat2d_ref = ref2d.as_numpy_array()
print dat2d_ref
from matplotlib import pyplot as plt

plt.imshow(dat2d_ref , interpolation = "nearest")
plt.show()


