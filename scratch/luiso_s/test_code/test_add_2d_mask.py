from dials.algorithms.integration import mask_add_2d
from scitbx.array_family import flex

nrow = ncol = 25

mask_01 = flex.int(flex.grid(nrow, ncol), 3)
mask_02 = flex.int(flex.grid(nrow, ncol), 3)


for row in range(int(nrow / 5), int(3 * nrow / 5)):
  for col in range(int(ncol / 5), int(3 * ncol / 5)):
    mask_01[row, col] = 5

for row in range(int(2 * nrow / 5), int(4 * nrow / 5)):
  for col in range(int(2 * ncol / 5), int(4 * ncol / 5)):
    mask_02[row, col] = 5

for diag in range(nrow):
  mask_01[diag, diag] = 0
  mask_02[diag, 13] = 0

from matplotlib import pyplot as plt

plt.imshow(mask_01.as_numpy_array(), interpolation = "nearest")
plt.show()
plt.imshow(mask_02.as_numpy_array(), interpolation = "nearest")
plt.show()


#sumation = mask_add_2d( mask_01, mask_02)
mask_01 = mask_add_2d( mask_01, mask_02)

plt.imshow(mask_01.as_numpy_array(), interpolation = "nearest")
plt.show()
plt.imshow(mask_02.as_numpy_array(), interpolation = "nearest")
plt.show()

from dials.scratch.luiso_s import write_2d_mask
write_2d_mask(mask_01)

#plt.imshow(sumation.as_numpy_array(), interpolation = "nearest")
#plt.show()
