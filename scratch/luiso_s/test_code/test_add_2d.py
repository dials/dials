from dials.algorithms.integration import add_2d
from scitbx.array_family import flex
from dials.scratch.luiso_s import model_2d

nrow = ncol = 25

ref2d_01 = model_2d(nrow, ncol, 3, 6, 0.2, 85, 0.5)
ref2d_02 = model_2d(nrow+1, ncol+1, 3, 6, 0.8, 85, 0.5)
descr = flex.double(flex.grid(1, 3))
descr[0, 0] = 12.5
descr[0, 1] = 12.5
descr[0, 2] = 2.5
#tmp_ref01 = ref2d_01[:,:]
#tmp_ref02 = ref2d_02[:,:]
print id(ref2d_01)
print id(ref2d_02)
print
sumation = add_2d(descr, ref2d_01, ref2d_02)
#sumation = add_2d(descr, tmp_ref01, tmp_ref02)
print
print id(ref2d_01)
print id(ref2d_02)
print id(sumation)


from matplotlib import pyplot as plt

plt.imshow(ref2d_01.as_numpy_array(), interpolation = "nearest")
plt.show()
plt.imshow(ref2d_02.as_numpy_array(), interpolation = "nearest")
plt.show()

plt.imshow(sumation.as_numpy_array(), interpolation = "nearest")
plt.show()
