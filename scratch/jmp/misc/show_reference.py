
from dials.model.serialize import load

reference = load.reference('/home/upc86896/Data/TRP_M1S3_2_/reference.pickle')

for i in range(len(reference)):
  from matplotlib import pylab
  profile = reference.profile(i)
  from scitbx.array_family import flex
  vmin, vmax = flex.min(profile), flex.max(profile)
  for j in range(1, 10):
    pylab.subplot(3, 3, j)
    pylab.imshow(profile.as_numpy_array()[j-1], vmin=vmin, vmax=vmax)
  pylab.show()
