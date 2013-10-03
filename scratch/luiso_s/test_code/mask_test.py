from __future__ import division
from dials.algorithms.shoebox import build_mask

x = build_mask()
#import numpy
print "x =", x.as_numpy_array()
