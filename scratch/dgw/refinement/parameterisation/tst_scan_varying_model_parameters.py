#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

from scan_varying_model_parameters import *
from math import pi

myparam = ScanVaryingParameterSet(1.0, 5)

#print dir (myparam)

# adjust one of the values
vals = myparam.value
vals[3] = 2.0
vals[4] = 2.0
myparam.value = vals
print myparam.value

# treat x range as the phi range of a scan
smoother = GaussianSmoother((0, 2 * pi), 3)

print smoother.num_values()
print smoother.num_samples()
print smoother.num_average()
print smoother.positions()
print smoother.spacing()

smooth_vals = [2*pi*e/100 for e in range(1, 100)]
for e in smooth_vals:

    val, weights, sumweights = smoother.value_weight(e, myparam)
    print e, val,
    for i in weights: print i,
    print sumweights
#print pi, smoother.value_weight(pi, myparam)
