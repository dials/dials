from scitbx.array_family import flex
from dials.algorithms.peak_finding import ref_2d



flex_ref2d_g = ref_2d(25, 25, 5, 3, .33333, 25, 0.9)
flex_ref2d_l = ref_2d(25, 25, 5, 3, .33333, 25, 0.1)


dat2d_ref_g = flex_ref2d_g.as_numpy_array()
dat2d_ref_l = flex_ref2d_l.as_numpy_array()

print dat2d_ref_g
print dat2d_ref_l


from matplotlib import pyplot as plt

plt.imshow(dat2d_ref_g , interpolation = "nearest")
plt.show()

plt.imshow(dat2d_ref_l , interpolation = "nearest")
plt.show()


