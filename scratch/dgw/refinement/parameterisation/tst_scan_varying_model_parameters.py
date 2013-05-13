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

###########################################################
# Start by testing a bare parameter set with the smoother #
###########################################################

# 7 values, all set to 1.0
myparam = ScanVaryingParameterSet(1.0, 7)

# Adjust a couple of the values
myparam.value[3:4] = [2.0, 2.0]

# Make a smoother with x_range as an 'image range', between 1 and 100. This
# smoother needs 5 intervals (for 7 total values). The default smoother uses
# an averaging window of 3 values
smoother = GaussianSmoother((1, 100), 5)

assert smoother.num_values() == 7
assert smoother.num_samples() == 5
assert smoother.num_average() == 3

# The smoother positions depend on the number of intervals but not the x_range
assert smoother.positions() == [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5]

# By contrast the spacing is in the original, unnormalised coordinate
assert smoother.spacing() == 19.8

# FIXME remove this printout, replace with optional scatterplot
smooth_vals = [e for e in range(1, 101)]
for e in smooth_vals:

    val, weights, sumweights = smoother.value_weight(e, myparam)
    print e, val,
    for i in weights: print i,
    print sumweights

print "OK"

#################################################################
# Now test a full ScanVaryingCrystalOrientationParameterisation #
#################################################################

# Use a class to wrap up ScanVaryingCrystalOrientationParameterisation with
# get_state overridden, so it can be passed to the existing FD derivative code.
class TestModel(ScanVaryingCrystalOrientationParameterisation):

    def __init__(self, image_number, *args):
        self.image_number = image_number
        ScanVaryingCrystalOrientationParameterisation.__init__(self, *args)

    def set_time_point(self, t):
        self.image_number = t

    def get_state(self):

        '''override get state to do so only at the requested t'''

        return ScanVaryingCrystalOrientationParameterisation.get_state(self,
            self.image_number)


from dials.algorithms.refinement \
    import get_fd_gradients, dR_from_axis_and_angle, random_param_shift
import random
from libtbx.test_utils import approx_equal
from cctbx.uctbx import unit_cell
from scitbx import matrix
from dials.model.experiment.crystal_model import Crystal

def random_direction_close_to(vector):
    return vector.rotate_around_origin(matrix.col(
                (random.random(),
                 random.random(),
                 random.random())).normalize(),
                 random.gauss(0, 1.0),  deg = True)

# Make a random crystal
a = random.uniform(10,50) * random_direction_close_to(matrix.col((1, 0, 0)))
b = random.uniform(10,50) * random_direction_close_to(matrix.col((0, 1, 0)))
c = random.uniform(10,50) * random_direction_close_to(matrix.col((0, 0, 1)))
xl = Crystal(a, b, c)

# Let's say we have a scan of 100 images
image_range = (1, 100)

# Parameterise the crystal with the image range and five intervals. Use
# TestModel to explore gradients at image 50.
xl_op = TestModel(50, xl, image_range, 5)

# FIXME remove this printout, replace with optional scatterplot
print "testing smoothed values for phi1 parameter"
smooth_vals = [e for e in range(1, 101)]
for e in smooth_vals:

    val, weights, sumweights = smoother.value_weight(e, xl_op._param_sets[0])
    print e, val,
    for i in weights: print i,
    print sumweights
print

# How many parameters?
num_param = xl_op.num_free()

# shift the parameters away from zero
p_vals = xl_op.get_p()
print "Original parameter vals", p_vals
sigmas = [0.1] * len(p_vals)
new_vals = random_param_shift(p_vals, sigmas)
xl_op.set_p(new_vals)
p_vals = xl_op.get_p()
print "Shifted parameter vals", p_vals

# compare analytical and finite difference derivatives at image 50
an_ds_dp = xl_op.get_ds_dp(50)
fd_ds_dp = get_fd_gradients(xl_op, [1.e-6 * pi/180] * num_param)
pnames = xl_op.get_pnames()

null_mat = matrix.sqr((0., 0., 0., 0., 0., 0., 0., 0., 0.))
for e, f in zip(an_ds_dp, fd_ds_dp):
    assert(approx_equal((e - f), null_mat, eps = 1.e-6))

print "OK"

#################################
# Test a few random time points #
#################################
for t in [random.uniform(1, 100) for x in range(50)]:

    xl_op.set_time_point(t)
    an_ds_dp = xl_op.get_ds_dp(t)
    fd_ds_dp = get_fd_gradients(xl_op, [1.e-6 * pi/180] * num_param)
    #print t
    #print "Gradients:"
    #for s, a, f in zip(pnames, an_ds_dp, fd_ds_dp):
    #    print s
    #    print a
    #    print f
    #    print "diff:", a-f
    #    print
    #
    for e, f in zip(an_ds_dp, fd_ds_dp):
        assert(approx_equal((e - f), null_mat, eps = 1.e-6))

print "OK"

######################################################################
# Random initial orientations, random parameter shifts, random times #
######################################################################

attempts = 100
failures = 0
for i in range(attempts):

    # make a random crystal and parameterise it
    a = random.uniform(10,50) * random_direction_close_to(matrix.col((1, 0, 0)))
    b = random.uniform(10,50) * random_direction_close_to(matrix.col((0, 1, 0)))
    c = random.uniform(10,50) * random_direction_close_to(matrix.col((0, 0, 1)))
    xl = Crystal(a, b, c)
    xl_op = TestModel(50, xl, image_range, 5)

    # apply random parameter shifts to the orientation (2.0 mrad each checkpoint)
    p_vals = xl_op.get_p()
    sigmas = [2.0] * len(p_vals)
    new_vals = random_param_shift(p_vals, sigmas)
    xl_op.set_p(new_vals)

    # select random time point at which to make comparisons
    t = random.uniform(*image_range)
    xl_op.set_time_point(t)

    # compare analytical and finite difference derivatives
    xl_op_an_ds_dp = xl_op.get_ds_dp(t)
    xl_op_fd_ds_dp = get_fd_gradients(xl_op, [1.e-6 * pi/180] * num_param)

    for j in range(num_param):
        try:
            assert(approx_equal((xl_op_fd_ds_dp[j] - xl_op_an_ds_dp[j]),
                                null_mat, eps = 1.e-6))
        except Exception:
            failures += 1
            print "for try", i
            print "failure for parameter number", j
            print "of the orientation parameterisation"
            print "with fd_ds_dp = "
            print fd_ds_dp[j]
            print "and an_ds_dp = "
            print an_ds_dp[j]
            print "so that difference fd_ds_dp - an_ds_dp ="
            print fd_ds_dp[j] - an_ds_dp[j]

if failures == 0: print "OK"
