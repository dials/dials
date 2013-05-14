#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

# Python imports
from __future__ import division
from math import pi
import random

# CCTBX imports
from libtbx.test_utils import approx_equal
from scitbx import matrix

# DIALS imports
from dials.algorithms.refinement \
    import get_fd_gradients, random_param_shift
from scan_varying_model_parameters import *
from dials.model.experiment.crystal_model import Crystal

class SmootherTest(object):

    '''Test a bare parameter set with the smoother'''
    def __init__(self, plots = False):

        # make scatterplots
        self.do_plots = plots

        # 7 values, all set to 1.0
        self.myparam = ScanVaryingParameterSet(1.0, 7)

        # Adjust a couple of the values
        self.myparam.value[3:4] = [2.0, 2.0]

        # Make a smoother with x_range as an 'image range', between 1 and 100.
        # This smoother needs 5 intervals (for 7 total values). The default
        # smoother uses an averaging window of 3 values
        self.smoother = GaussianSmoother((1, 100), 5)

    def run(self):

        assert self.smoother.num_values() == 7
        assert self.smoother.num_samples() == 5
        assert self.smoother.num_average() == 3

        # The smoother positions depend on the number of intervals but not on
        # the x_range
        assert self.smoother.positions() == [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5]

        # By contrast the spacing is in units of the original, unnormalised
        # coordinate
        assert self.smoother.spacing() == 19.8

        smooth_at = [e for e in range(1, 101)]
        data = [self.smoother.value_weight(e, self.myparam) for e in smooth_at]
        vals = [v for v, w, sw in data]
        assert len(smooth_at) == len(vals)

        if self.do_plots:
            try:
                import matplotlib.pyplot as plt
                plt.ion()
                plt.scatter(smooth_at, vals)
                plt.draw()
            except ImportError as e:
                print "pyplot not available", e

        print "OK"


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


class TestScanVaryingCrystalOrientationParameterisation(object):

    '''Test a ScanVaryingCrystalOrientationParameterisation'''

    def __init__(self, plots = False):

        # Do we make scatterplots?
        self.do_plots = plots

        # Let's say we have a scan of 100 images
        self.image_range = (1, 100)

        # Make a random crystal
        a = random.uniform(10,50) * \
            self.random_direction_close_to(matrix.col((1, 0, 0)))
        b = random.uniform(10,50) * \
            self.random_direction_close_to(matrix.col((0, 1, 0)))
        c = random.uniform(10,50) * \
            self.random_direction_close_to(matrix.col((0, 0, 1)))
        self.xl = Crystal(a, b, c)

    def random_direction_close_to(self, vector):
        return vector.rotate_around_origin(matrix.col(
                    (random.random(),
                     random.random(),
                     random.random())).normalize(),
                     random.gauss(0, 1.0),  deg = True)

    def test_num_intervals(self, nintervals):
        '''Test a range of different numbers of intervals'''

        # Parameterise the crystal with the image range and five intervals. Init
        # TestModel to explore gradients at image 50, but actually we will try
        # various time points.
        xl_op = TestModel(50, self.xl, self.image_range, nintervals)

        # How many parameters?
        num_param = xl_op.num_free()

        # shift the parameters away from zero
        p_vals = xl_op.get_p()
        sigmas = [0.1] * len(p_vals)
        new_vals = random_param_shift(p_vals, sigmas)
        xl_op.set_p(new_vals)
        p_vals = xl_op.get_p()
        #print "Shifted parameter vals", p_vals

        # compare analytical and finite difference derivatives at image 50
        an_ds_dp = xl_op.get_ds_dp(50)
        fd_ds_dp = get_fd_gradients(xl_op, [1.e-6 * pi/180] * num_param)
        pnames = xl_op.get_pnames()

        null_mat = matrix.sqr((0., 0., 0., 0., 0., 0., 0., 0., 0.))
        for e, f in zip(an_ds_dp, fd_ds_dp):
            assert(approx_equal((e - f), null_mat, eps = 1.e-6))

        # Now test gradients at equally spaced time points across the whole
        # range
        num_points = 50
        smooth_at = []
        phi1_data = []
        phi2_data = []
        phi3_data = []
        step_size = (self.image_range[1] - self.image_range[0]) / num_points
        for t in [self.image_range[0] + e * step_size \
                    for e in range(num_points + 1)]:

            # collect data for plot
            smooth_at.append(t)
            phi1_data.append(xl_op._smoother.value_weight(t, xl_op._param_sets[0])[0])
            phi2_data.append(xl_op._smoother.value_weight(t, xl_op._param_sets[1])[0])
            phi3_data.append(xl_op._smoother.value_weight(t, xl_op._param_sets[2])[0])

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

        if self.do_plots:
            try:
                import matplotlib.pyplot as plt
                plt.ion()
                plt.clf()
                plt.subplot(311)
                plt.cla()
                plt.scatter(smooth_at, phi1_data)
                plt.title("Phi1")
                plt.xlabel("image number")
                plt.ylabel("Phi1 (mrad)")
                plt.subplot(312)
                plt.cla()
                plt.scatter(smooth_at, phi2_data)
                plt.title("Phi2")
                plt.xlabel("image number")
                plt.ylabel("Phi2 (mrad)")
                plt.subplot(313)
                plt.cla()
                plt.scatter(smooth_at, phi3_data)
                plt.title("Phi3")
                plt.xlabel("image number")
                plt.ylabel("Phi3 (mrad)")
                plt.suptitle("Parameter smoothing with %d intervals" % nintervals)
                plt.draw()
            except ImportError as e:
                print "pyplot not available", e

        print "OK"

    def test_random(self):

        '''Test random initial orientations, random parameter shifts and random
        times'''

        attempts = 100
        failures = 0
        null_mat = matrix.sqr((0., 0., 0., 0., 0., 0., 0., 0., 0.))

        for i in range(attempts):

            # make a new random crystal and parameterise it
            a = random.uniform(10,50) * \
                    self.random_direction_close_to(matrix.col((1, 0, 0)))
            b = random.uniform(10,50) * \
                    self.random_direction_close_to(matrix.col((0, 1, 0)))
            c = random.uniform(10,50) * \
                    self.random_direction_close_to(matrix.col((0, 0, 1)))
            xl = Crystal(a, b, c)

            xl_op = TestModel(50, xl, self.image_range, 5)

            # How many parameters?
            num_param = xl_op.num_free()

            # apply random parameter shifts to the orientation (2.0 mrad each
            # checkpoint)
            p_vals = xl_op.get_p()
            sigmas = [2.0] * len(p_vals)
            new_vals = random_param_shift(p_vals, sigmas)
            xl_op.set_p(new_vals)

            # select random time point at which to make comparisons
            t = random.uniform(*self.image_range)
            xl_op.set_time_point(t)

            # compare analytical and finite difference derivatives
            xl_op_an_ds_dp = xl_op.get_ds_dp(t)
            xl_op_fd_ds_dp = get_fd_gradients(xl_op, [1.e-6 * pi/180] * \
                                              num_param)

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
                    print xl_op_fd_ds_dp[j]
                    print "and an_ds_dp = "
                    print xl_op_an_ds_dp[j]
                    print "so that difference fd_ds_dp - an_ds_dp ="
                    print xl_op_fd_ds_dp[j] - xl_op_an_ds_dp[j]

        if failures == 0: print "OK"

    def run(self):

        for n in (1, 2, 3, 4, 5, 6, 7):
            self.test_num_intervals(n)

        self.test_random()


if __name__ == '__main__':

    test = SmootherTest(plots = False)
    test.run()

    test = TestScanVaryingCrystalOrientationParameterisation(plots = False)
    test.run()
