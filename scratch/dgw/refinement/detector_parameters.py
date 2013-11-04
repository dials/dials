#!/usr/bin/env python

# Copyright (C) (2012) David Waterman, STFC Rutherford Appleton Laboratory, UK.
# This code is developed as part of the DIALS project and is provided for
# testing purposes only

from __future__ import division
from model_parameters import Parameter, ModelParameterisation
from scitbx import matrix
from dials.model.experiment import Panel, Detector
from dials.model.experiment import detector_factory
from math import sin, cos, pi, sqrt
from dials.scratch.dgw.refinement \
    import dR_from_axis_and_angle, get_fd_gradients, random_param_shift

class DetectorParameterisationSinglePanel(ModelParameterisation):
    """implementation of parameterisation for a single abstract panel
    plane, with angles expressed in mrad"""

    def __init__(self, detector):

        # the state of the detector model is comprised of:
        #
        # * a vector 'dorg' locating the origin of the detector model
        # * two orthogonal unit directions 'd1' and 'd2' forming a plane
        #   with its origin at the end of the vector dorg
        # * a third unit direction 'dn', orthogonal to both 'd1' & 'd2'.
        # * offsets to locate the origin of each sensor frame from the
        #   tip of the dorg vector, in terms of the coordinate system
        #   formed by d1, d2 and dn.
        #
        # Held separately in attribute 'models' are:
        # * references to the sensor objects contained in this model
        #
        # For this simplified class there is only a single sensor frame
        # and the vector dn is not required, because the plane formed by
        # d1 and d2 is coplanar with the sensor plane. Therefore the
        # offset is fully in terms of d1 and d2

        # set up the initial state of the detector model from the
        # orientation of the single sensor it contains

        # get some vectors we need from the Panel
        panel = detector[0]
        so = matrix.col(panel.get_origin())
        d1 = matrix.col(panel.get_fast_axis())
        d2 = matrix.col(panel.get_slow_axis())
        dn = matrix.col(panel.get_normal())

        # we choose the dorg vector to terminate in the centre of the
        # sensor, and the offset between the end of the dorg vector and
        # the sensor origin is a coordinate matrix with elements in the
        # basis d1, d2, dn
        panel_lim = panel.get_image_size_mm()
        offset = matrix.col((-1. * panel_lim[0] / 2.,
                             -1. * panel_lim[1] / 2.,
                              0.))

        dorg = so - offset[0] * d1 - offset[1] * d2

        # Set up the initial state. There are multiple items of interest, so
        # use a dictionary here (note for a single sensor model we can do
        # without the first 3 of these, but will need them for multiple sensors)
        istate = {'d1':d1,
                  'd2':d2,
                  'dn':dn,
                  'offset':offset}

        # set up the parameters.
        # distance from lab origin to detector model plane along its
        # normal, in initial orientation
        distance = panel.get_distance()
        dist = Parameter(distance, dn, 'length')

        # shift in the detector model plane to locate dorg, in initial
        # orientation
        shift = dorg - dn * distance
        shift1 = Parameter(shift.dot(d1), d1, 'length')
        shift2 = Parameter(shift.dot(d2), d2, 'length')

        # rotations of the plane through its origin about:
        # 1) axis normal to initial orientation
        # 2) d1 axis of initial orientation
        # 3) d2 axis of initial orientation
        tau1 = Parameter(0, dn, 'angle')
        tau2 = Parameter(0, d1, 'angle')
        tau3 = Parameter(0, d2, 'angle')

        # build the parameter list in a specific,  maintained order
        p_list = [dist, shift1, shift2, tau1, tau2, tau3]

        # set up the list of model objects being parameterised (here
        # just the detector containing a single panel)
        models = [detector]

        # set up the base class
        ModelParameterisation.__init__(self, models, istate, p_list)

        # call compose to calculate all the derivatives
        self.compose()

    def compose(self):

        # extract items from the initial state
        id1 = self._initial_state['d1']
        id2 = self._initial_state['d2']
        ioffset = self._initial_state['offset']

        # extract parameters from the internal list
        dist, shift1, shift2, tau1, tau2, tau3 = self._param

        # convert angles to radians
        tau1rad = tau1.value / 1000.
        tau2rad = tau2.value / 1000.
        tau3rad = tau3.value / 1000.

        # compose rotation matrices and their first order derivatives
        Tau1 = (tau1.axis).axis_and_angle_as_r3_rotation_matrix(tau1rad,
                                                                deg=False)
        dTau1_dtau1 = dR_from_axis_and_angle(tau1.axis, tau1rad, deg=False)

        Tau2 = (tau2.axis).axis_and_angle_as_r3_rotation_matrix(tau2rad,
                                                                deg=False)
        dTau2_dtau2 = dR_from_axis_and_angle(tau2.axis, tau2rad, deg=False)

        Tau3 = (tau3.axis).axis_and_angle_as_r3_rotation_matrix(tau3rad,
                                                                deg=False)
        dTau3_dtau3 = dR_from_axis_and_angle(tau3.axis, tau3rad, deg=False)

        Tau32 = Tau3 * Tau2
        Tau321 = Tau32 * Tau1

        # Compose new state
        # =================
        # First the frame positioned at a distance from the lab origin
        P0 = dist.value * dist.axis # distance along initial detector normal
        Px = P0 + id1 # point at the end of d1 in lab frame
        Py = P0 + id2 # point at the end of d2 in lab frame

        # detector shift vector
        dsv = P0 + shift1.value * shift1.axis + shift2.value * shift2.axis

        # compose dorg point
        dorg = Tau321 * dsv - Tau32 * P0 + P0

        # compose d1, d2 and dn and ensure frame remains orthonormal.
        d1 = (Tau321 * (Px - P0)).normalize()
        d2 = (Tau321 * (Py - P0)).normalize()
        dn = d1.cross(d2).normalize()
        # NB dn not actually used in this simple model; calculation left
        # here as a reminder for future extension
        d2 = dn.cross(d1)

        # compose new sensor origin
        o = dorg + ioffset[0] * d1 + ioffset[1] * d2

        # compose new sensor dir1, dir2. For this model with a single
        # sensor, this is trivial as they are the same as d1 and d2
        dir1 = d1
        dir2 = d2

        # now update the panel with its new position and orientation.
        # The detector is the first model in _models, the panel is the
        # first in the detector
        (self._models[0])[0].set_local_frame(dir1, dir2, o)

        # calculate derivatives of the state wrt parameters
        # =================================================
        # Start with the dorg vector, where
        # dorg = Tau321 * dsv - Tau32 * P0 + P0

        # derivative wrt dist
        dP0_ddist = dist.axis
        ddsv_ddist = dP0_ddist
        ddorg_ddist = Tau321 * ddsv_ddist - Tau32 * dP0_ddist + \
                      dP0_ddist

        # derivative wrt shift1
        ddsv_dshift1 = shift1.axis
        ddorg_dshift1 = Tau321 * ddsv_dshift1

        # derivative wrt shift2
        ddsv_dshift2 = shift2.axis
        ddorg_dshift2 = Tau321 * ddsv_dshift2

        # derivative wrt tau1
        dTau321_dtau1 = Tau32 * dTau1_dtau1
        ddorg_dtau1 = dTau321_dtau1 * dsv

        # derivative wrt tau2
        dTau32_dtau2 = Tau3 * dTau2_dtau2
        dTau321_dtau2 = dTau32_dtau2 * Tau1
        ddorg_dtau2 = dTau321_dtau2 * dsv - dTau32_dtau2 * P0

        # derivative wrt tau3
        dTau32_dtau3 = dTau3_dtau3 * Tau2
        dTau321_dtau3 = dTau32_dtau3 * Tau1
        ddorg_dtau3 = dTau321_dtau3 * dsv - dTau32_dtau3 * P0

        # Now derivatives of the direction d1, where
        # d1 = (Tau321 * (Px - P0)).normalize()
        # For calc of derivatives ignore the normalize(), which should
        # be unnecessary anyway as Px - P0 is a unit vector and Tau321 a
        # pure rotation.

        # derivative wrt dist
        # dPx_ddist = dist.axis; dP0_ddist = dist.axis, so these cancel
        dd1_ddist = matrix.col((0., 0., 0.))

        # derivative wrt shift1
        dd1_dshift1 = matrix.col((0., 0., 0.))

        # derivative wrt shift2
        dd1_dshift2 = matrix.col((0., 0., 0.))

        # derivative wrt tau1
        dd1_dtau1 = dTau321_dtau1 * (Px - P0)

        # derivative wrt tau2
        dd1_dtau2 = dTau321_dtau2 * (Px - P0)

        # derivative wrt tau3
        dd1_dtau3 = dTau321_dtau3 * (Px - P0)

        # Derivatives of the direction d2, where
        # d2 = (Tau321 * (Py - P0)).normalize()

        # derivative wrt dist
        dd2_ddist = matrix.col((0., 0., 0.))

        # derivative wrt shift1
        dd2_dshift1 = matrix.col((0., 0., 0.))

        # derivative wrt shift2
        dd2_dshift2 = matrix.col((0., 0., 0.))

        # derivative wrt tau1
        dd2_dtau1 = dTau321_dtau1 * (Py - P0)

        # derivative wrt tau2
        dd2_dtau2 = dTau321_dtau2 * (Py - P0)

        # derivative wrt tau3
        dd2_dtau3 = dTau321_dtau3 * (Py - P0)

        # Derivatives of the direction dn, where
        # dn = d1.cross(d2).normalize()

        # derivative wrt dist
        ddn_ddist = matrix.col((0., 0., 0.))

        # derivative wrt shift1
        ddn_dshift1 = matrix.col((0., 0., 0.))

        # derivative wrt shift2
        ddn_dshift2 = matrix.col((0., 0., 0.))

        # derivative wrt tau1. Product rule for cross product applies
        ddn_dtau1 = dd1_dtau1.cross(d2) + dd2_dtau1.cross(d1)

        # derivative wrt tau2
        ddn_dtau2 = dd1_dtau2.cross(d2) + dd2_dtau2.cross(d1)

        # derivative wrt tau3
        ddn_dtau3 = dd1_dtau3.cross(d2) + dd2_dtau3.cross(d1)

        # calculate derivatives of the attached sensor matrix
        # ===================================================
        # sensor origin:
        # o = dorg + ioffset[0] * d1 + ioffset[1] * d2

        # derivative wrt dist
        do_ddist = ddorg_ddist + ioffset[0] * dd1_ddist + \
                   ioffset[1] * dd2_ddist

        # derivative wrt shift1
        do_dshift1 = ddorg_dshift1 + ioffset[0] * dd1_dshift1 + \
                     ioffset[1] * dd2_dshift1

        # derivative wrt shift2
        do_dshift2 = ddorg_dshift2 + ioffset[0] * dd1_dshift2 + \
                     ioffset[1] * dd2_dshift2

        # derivative wrt tau1
        do_dtau1 = ddorg_dtau1 + ioffset[0] * dd1_dtau1 + \
                   ioffset[1] * dd2_dtau1

        # derivative wrt tau2
        do_dtau2 = ddorg_dtau2 + ioffset[0] * dd1_dtau2 + \
                   ioffset[1] * dd2_dtau2

        # derivative wrt tau3
        do_dtau3 = ddorg_dtau3 + ioffset[0] * dd1_dtau3 + \
                   ioffset[1] * dd2_dtau3

        # sensor dir1:
        # dir1 = d1

        # derivative wrt dist
        ddir1_ddist = dd1_ddist

        # derivative wrt shift1
        ddir1_dshift1 = dd1_dshift1

        # derivative wrt shift2
        ddir1_dshift2 = dd1_dshift2

        # derivative wrt tau1
        ddir1_dtau1 = dd1_dtau1

        # derivative wrt tau2
        ddir1_dtau2 = dd1_dtau2

        # derivative wrt tau3
        ddir1_dtau3 = dd1_dtau3

        # sensor dir2:
        # dir2 = d2

        # derivative wrt dist
        ddir2_ddist = dd2_ddist

        # derivative wrt shift1
        ddir2_dshift1 = dd2_dshift1

        # derivative wrt shift2
        ddir2_dshift2 = dd2_dshift2

        # derivative wrt tau1
        ddir2_dtau1 = dd2_dtau1

        # derivative wrt tau2
        ddir2_dtau2 = dd2_dtau2

        # derivative wrt tau3
        ddir2_dtau3 = dd2_dtau3

        # combine these vectors together into derivatives of the sensor
        # matrix d and store them, converting angles back to mrad

        # derivative wrt dist
        self._dstate_dp[0] = matrix.sqr(ddir1_ddist.elems +
                              ddir2_ddist.elems +
                              do_ddist.elems).transpose()

        # derivative wrt shift1
        self._dstate_dp[1] = matrix.sqr(ddir1_dshift1.elems +
                              ddir2_dshift1.elems +
                              do_dshift1.elems).transpose()

        # derivative wrt shift2
        self._dstate_dp[2] = matrix.sqr(ddir1_dshift2.elems +
                              ddir2_dshift2.elems +
                              do_dshift2.elems).transpose()

        # derivative wrt tau1
        self._dstate_dp[3] = matrix.sqr(ddir1_dtau1.elems +
                              ddir2_dtau1.elems +
                              do_dtau1.elems).transpose() / 1000.

        # derivative wrt tau2
        self._dstate_dp[4] = matrix.sqr(ddir1_dtau2.elems +
                              ddir2_dtau2.elems +
                              do_dtau2.elems).transpose() / 1000.

        # derivative wrt tau3
        self._dstate_dp[5] = matrix.sqr(ddir1_dtau3.elems +
                              ddir2_dtau3.elems +
                              do_dtau3.elems).transpose() / 1000.

        return

    def get_state(self):
        # only a single panel
        panel = (self._models[0])[0]
        return matrix.sqr(panel.get_d_matrix())

def random_panel(lim = (0, 50)):
    """For testing, return a square panel with a randomised position
    and orientation"""

    import random

    # start with a randomised origin vector
    o = matrix.col((random.uniform(-200, 200),
                    random.uniform(-200, 200),
                    random.uniform(-200, 200)))

    # two orthogonal unit vectors randomly oriented in the normal plane
    # of the origin vector
    u1 = o.ortho().normalize()
    u2 = o.cross(u1).normalize()
    theta = random.uniform(0, 2. * pi)
    u1 = u1.rotate_around_origin(o, pi/12)
    u2 = u2.rotate_around_origin(o, pi/12)

    # offset the plane normal from the origin vector by random rotations
    # of up to 45 degrees for each direction
    u1 = u1.rotate_around_origin(u2, random.uniform(-pi/2., pi/2.))
    u2 = u2.rotate_around_origin(u1, random.uniform(-pi/2., pi/2.))

    return Panel("PAD", u1, u2, o,
            (lim[1]/200, lim[1]/200), (200, 200), (0, 2e20))

if __name__ == '__main__':

    from libtbx.test_utils import approx_equal

    # set up a simple detector frame with directions aligned with
    # principal axes and sensor origin located on the z-axis at -110
    d1 = matrix.col((1, 0, 0))
    d2 = matrix.col((0, -1, 0))
    #lim = (0,50)
    npx_fast = 1475
    npx_slow = 1679
    pix_size_f = pix_size_s = 0.172
    detector = detector_factory.make_detector("PAD", d1, d2,
        matrix.col((0, 0, -110)), (pix_size_f, pix_size_s),
        (npx_fast, npx_slow), (0, 2e20))

    dp = DetectorParameterisationSinglePanel(detector)

    # Test change of parameters
    # =========================

    # 1. shift detector plane so that the z-axis intercepts its centre
    # at a distance of 100 along the initial normal direction. As the
    # initial normal is along -z, we expect the frame to intercept the
    # z-axis at -100.

    p_vals = dp.get_param_vals()
    p_vals[0:3] = [100., 0., 0.]
    dp.set_param_vals(p_vals)
    panel = dp._models[0]
    v1 = matrix.col(panel.get_origin())
    v2 = matrix.col((0., 0., 1.))
    assert(approx_equal(v1.dot(v2), -100.))

    # 2. rotate frame around its initial normal by +90 degrees. Only d1
    # and d2 should change. As we rotate clockwise around the initial
    # normal (-z direction) then d1 should rotate onto the original
    # direction d2, and d2 should rotate to negative of the original
    # direction d1

    p_vals[3] = 1000. * pi/2 # set tau1 value
    dp.set_param_vals(p_vals)

    panel = dp._models[0]
    assert(approx_equal(matrix.col(panel.get_fast_axis()).dot(dp._initial_state['d1']), 0.))
    assert(approx_equal(matrix.col(panel.get_slow_axis()).dot(dp._initial_state['d2']), 0.))
    assert(approx_equal(matrix.col(panel.get_normal()).dot(dp._initial_state['dn']), 1.))

    # 3. no rotation around initial normal, +10 degrees around initial
    # d1 direction and +10 degrees around initial d2. Check d1 and d2
    # match paper calculation

    p_vals[3] = 0.    # tau1
    p_vals[4] = 1000. * pi/18 # tau2
    p_vals[5] = 1000. * pi/18 # tau3
    dp.set_param_vals(p_vals)

    # Debugging - left here commented out as an example
    #print "initial state"
    #print "d1 = (%.3f, %.3f, %.3f)" % dp._initial_state['d1'].elems
    #print "d2 = (%.3f, %.3f, %.3f)" % dp._initial_state['d2'].elems
    #print "dn = (%.3f, %.3f, %.3f)" % dp._initial_state['dn'].elems
    #print "composed state"
    #print "dir1 = (%.3f, %.3f, %.3f)" % dp._models[0].dir1
    #print "d2 = (%.9f, %.9f, %.9f)" % dp._models[0].dir2
    #print "dn = (%.3f, %.3f, %.3f)" % dp._models[0].normal

    # paper calculation values
    v1 = matrix.col((cos(pi/18), 0, sin(pi/18)))
    v2 = matrix.col((sin(pi/18)**2,
                     -cos(pi/18),
                     sqrt((2*sin(pi/36)*sin(pi/18))**2 - sin(pi/18)**4) - sin(pi/18)))

    panel = dp._models[0]
    assert(approx_equal(matrix.col(panel.get_fast_axis()).dot(v1), 1.))
    assert(approx_equal(matrix.col(panel.get_slow_axis()).dot(v2), 1.))

    # 4. Test fixing and unfixing of parameters
    p_vals = [100., 0., 0., 1000.*pi/18, 1000.*pi/18, 1000.*pi/18]
    dp.set_param_vals(p_vals)
    f = dp.get_fixed()
    f[0:3] = [True] * 3
    dp.set_fixed(f)
    p_vals2 = [0., 0., 0.]
    dp.set_param_vals(p_vals2)
    assert(dp.get_param_vals(only_free = False) == [100., 0., 0., 0., 0., 0.])

    an_ds_dp = dp.get_ds_dp()
    assert(len(an_ds_dp) == 3)

    f[0:3] = [False] * 3
    dp.set_fixed(f)
    p_vals = dp.get_param_vals()
    p_vals2 = [a + b for a, b in zip(p_vals, [-10., 1., 1., 0., 0., 0.])]
    dp.set_param_vals(p_vals2)
    assert(dp.get_param_vals() == [90., 1., 1., 0., 0., 0.])

    # 5. Tests of the calculation of derivatives

    # Now using parameterisation in mrad

    # random initial orientations with a random parameter shift at each
    attempts = 100
    failures = 0
    for i in range(attempts):

        # create random initial position
        dp = DetectorParameterisationSinglePanel(Detector(
            random_panel()))

        # apply a random parameter shift
        p_vals = random_param_shift(p_vals, [10, 10, 10, 1000.*pi/18,
                                             1000.*pi/18, 1000.*pi/18])
        dp.set_param_vals(p_vals)

        # obtain current sensor state
        #state = dp.get_state()

        # compare analytical and finite difference derivatives
        an_ds_dp = dp.get_ds_dp()
        fd_ds_dp = get_fd_gradients(dp, [1.e-6] * 3 + [1.e-4 * pi/180] * 3)

        for j in range(6):
            try:
                assert(approx_equal((fd_ds_dp[j] - an_ds_dp[j]),
                        matrix.sqr((0., 0., 0.,
                                    0., 0., 0.,
                                    0., 0., 0.)), eps = 1.e-6))
            except Exception:
                failures += 1
                print "for try", i
                print "failure for parameter number", j
                print "with fd_ds_dp = "
                print fd_ds_dp[j]
                print "and an_ds_dp = "
                print an_ds_dp[j]
                print "so that difference fd_ds_dp - an_ds_dp ="
                print fd_ds_dp[j] - an_ds_dp[j]

    if failures == 0: print "OK"
