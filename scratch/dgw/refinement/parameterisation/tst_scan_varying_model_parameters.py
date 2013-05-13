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

myparam = ScanVaryingParameterSet(1.0, 7)

#print dir (myparam)

# adjust one of the values
vals = myparam.value
vals[3] = 2.0
vals[4] = 2.0
myparam.value = vals
print myparam.value

# treat x range as image range, between 1 and 100
smoother = GaussianSmoother((1, 100), 5)

print smoother.num_values()
print smoother.num_samples()
print smoother.num_average()
print smoother.positions()
print smoother.spacing()

smooth_vals = [e for e in range(1, 101)]
for e in smooth_vals:

    val, weights, sumweights = smoother.value_weight(e, myparam)
    print e, val,
    for i in weights: print i,
    print sumweights

class TestModel(ScanVaryingCrystalOrientationParameterisation):

    def __init__(self, image_number, *args):
        self.image_number = image_number
        ScanVaryingCrystalOrientationParameterisation.__init__(self, *args)

    def set_time_point(self, t):
        self.image_number = t

    def get_state(self):

        '''override get state to do so only at the requested t'''

        return ScanVaryingCrystalOrientationParameterisation.get_state(self, self.image_number)



if __name__ == '__main__':

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

    # make a random crystal and parameterise it
    a = random.uniform(10,50) * random_direction_close_to(matrix.col((1, 0, 0)))
    b = random.uniform(10,50) * random_direction_close_to(matrix.col((0, 1, 0)))
    c = random.uniform(10,50) * random_direction_close_to(matrix.col((0, 0, 1)))
    xl = Crystal(a, b, c)

    # Let's say we have a scan of 100 images
    image_range = (1, 100)

    # Parameterise with the image range and five intervals. Use the test to
    # epxlore gradients at image 50
    #xl_op = ScanVaryingCrystalOrientationParameterisation(xl, image_range, 5)
    xl_op = TestModel(50, xl, image_range, 5)
    #xl_ucp = CrystalUnitCellParameterisation(xl)

    print "testing smoothed values for phi1 parameter"
    smooth_vals = [e for e in range(1, 101)]
    for e in smooth_vals:

        val, weights, sumweights = smoother.value_weight(e, xl_op._param_sets[0])
        print e, val,
        for i in weights: print i,
        print sumweights
    print

    # how many parameters?
    num_param = xl_op.num_free()

    null_mat = matrix.sqr((0., 0., 0., 0., 0., 0., 0., 0., 0.))

    # shift the parameters away from zero
    p_vals = xl_op.get_p()
    print "original parameter vals", p_vals
    sigmas = [0.1] * len(p_vals)
    new_vals = random_param_shift(p_vals, sigmas)
    xl_op.set_p(new_vals)
    p_vals = xl_op.get_p()
    print "shifted parameter vals", p_vals


    # compare analytical and finite difference derivatives at image 50
    an_ds_dp = xl_op.get_ds_dp(50)
    fd_ds_dp = get_fd_gradients(xl_op, [1.e-6 * pi/180] * num_param)
    pnames = xl_op.get_pnames()

    #print "Gradients:"
    #for s, a, f in zip(pnames, an_ds_dp, fd_ds_dp):
    #    print s
    #    print a
    #    print f
    #    print "diff:", a-f
    #    print

    for e, f in zip(an_ds_dp, fd_ds_dp):
        assert(approx_equal((e - f), null_mat, eps = 1.e-6))

    # try the same at various other time points
    for t in (1.0, 92, 100):

        xl_op.set_time_point(t)
        an_ds_dp = xl_op.get_ds_dp(t)
        fd_ds_dp = get_fd_gradients(xl_op, [1.e-6 * pi/180] * num_param)
        print t
        print "Gradients:"
        for s, a, f in zip(pnames, an_ds_dp, fd_ds_dp):
            print s
            print a
            print f
            print "diff:", a-f
            print

        for e, f in zip(an_ds_dp, fd_ds_dp):
            assert(approx_equal((e - f), null_mat, eps = 1.e-6))
    1/0

    #an_ds_dp = xl_ucp.get_ds_dp()
    #fd_ds_dp = get_fd_gradients(xl_ucp, [1.e-7] * xl_ucp.num_free())
    #for e, f in zip(an_ds_dp, fd_ds_dp):
    #    print e
    #    print f
    #    assert(approx_equal((e - f), null_mat, eps = 1.e-6))

    # random initial orientations with a random parameter shift at each
    attempts = 100
    failures = 0
    for i in range(attempts):

        # make a random crystal and parameterise it
        a = random.uniform(10,50) * random_direction_close_to(matrix.col((1, 0, 0)))
        b = random.uniform(10,50) * random_direction_close_to(matrix.col((0, 1, 0)))
        c = random.uniform(10,50) * random_direction_close_to(matrix.col((0, 0, 1)))
        xl = Crystal(a, b, c)
        xl_op = CrystalOrientationParameterisation(xl)
        xl_uc = CrystalUnitCellParameterisation(xl)

        # apply a random parameter shift to the orientation
        p_vals = xl_op.get_p()
        p_vals = random_param_shift(p_vals, [1000*pi/9, 1000*pi/9,
                                             1000*pi/9])
        xl_op.set_p(p_vals)

        # compare analytical and finite difference derivatives
        xl_op_an_ds_dp = xl_op.get_ds_dp()
        xl_op_fd_ds_dp = get_fd_gradients(xl_op, [1.e-5 * pi/180] * 3)

        # apply a random parameter shift to the unit cell

        print "\nCYCLE", i, "\n"
        print "apply random parameter shift"
        p_vals = xl_uc.get_p()
        cell_params = xl.get_unit_cell().parameters()
        print "old unit cell",cell_params
        cell_params = random_param_shift(cell_params, [1.] * 6)
        new_uc = unit_cell(cell_params)
        print "new unit cell",cell_params
        newB = matrix.sqr(new_uc.fractionalization_matrix()).transpose()
        S = symmetrize_reduce_enlarge(xl.get_space_group())
        S.set_orientation(orientation=newB)
        X = S.forward_independent_parameters()
        print "old p_vals", p_vals
        print "new p_vals", X, "\n"
        print "set these parameters in the model parameterisation"
        xl_uc.set_p(X)

        xl_uc_an_ds_dp = xl_ucp.get_ds_dp()
        print "\nnow doing finite differences about each parameter in turn"
        xl_uc_fd_ds_dp = get_fd_gradients(xl_ucp, [1.e-7] * xl_ucp.num_free())

        for j in range(3):
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

        for j in range(xl_ucp.num_free()):
            try:
                assert(approx_equal((xl_uc_fd_ds_dp[j] - xl_uc_an_ds_dp[j]),
                                    null_mat, eps = 1.e-6))
            except Exception:
                failures += 1
                print "for try", i
                print "failure for parameter number", j
                print "of the unit cell parameterisation"
                print "with fd_ds_dp = "
                print fd_ds_dp[j]
                print "and an_ds_dp = "
                print an_ds_dp[j]
                print "so that difference fd_ds_dp - an_ds_dp ="
                print fd_ds_dp[j] - an_ds_dp[j]

    if failures == 0: print "OK"
