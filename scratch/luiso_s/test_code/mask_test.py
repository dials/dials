from __future__ import division
from dials.algorithms.shoebox import build_mask

x = build_mask(23, 17, 3, 2, 8)
#import numpy
print "mask ="
print x.as_numpy_array()
