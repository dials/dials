"""Auxiliary functions for the refinement package"""

from __future__ import division
from math import sin, cos
from scitbx import matrix
from dials_refinement_ext import *
import random

def dR_from_axis_and_angle(axis, angle, deg=False):
    '''return the first derivative of a rotation matrix specified by its
    axis and angle'''

    # NB it is inefficient to do this separately from the calculation of
    # the rotation matrix itself, but it seems the Python interface to
    # scitbx does not have a suitable function. It might perhaps be
    # useful to write one, which could come straight from David Thomas'
    # RTMATS (present in Mosflm and MADNES).

    # NB RTMATS does calculation for a clockwise rotation of a vector
    # whereas axis_and_angle_as_r3_rotation_matrix does anticlockwise
    # rotation. Therefore flip the axis in here compared with
    # RTMATS in order to match the axis_and_angle_as_r3_rotation_matrix
    # convention

    assert axis.n in ((3,1), (1,3))
    if (deg): angle *= pi/180
    axis *= -1.
    ca, sa  = cos(angle), sin(angle)

    return(matrix.sqr((sa * axis[0] * axis[0] - sa ,
                    sa * axis[0] * axis[1] + ca * axis[2],
                    sa * axis[0] * axis[2] - ca * axis[1],
                    sa * axis[1] * axis[0] - ca * axis[2],
                    sa * axis[1] * axis[1] - sa,
                    sa * axis[1] * axis[2] + ca * axis[0],
                    sa * axis[2] * axis[0] + ca * axis[1],
                    sa * axis[2] * axis[1] - ca * axis[0],
                    sa * axis[2] * axis[2] - sa)))

def random_param_shift(vals, sigmas):
    '''Add a random (normal) shift to a parameter set, for testing'''

    assert len(vals) == len(sigmas)
    shifts = [random.gauss(0, sd) for sd in sigmas]
    newvals = [(x + y) for x, y in zip(vals, shifts)]

    return newvals

def get_fd_gradients(mp, deltas):
    '''Calculate centered finite difference gradients for each of the
    parameters of the model parameterisation mp.

    "deltas" must be a sequence of the same length as the parameter list, and
    contains the step size for the difference calculations for each parameter.
    '''

    #state = matrix.sqr(mp._models[0].d)

    p_vals = mp.get_p()
    assert len(deltas) == len(p_vals)
    fd_grad = []

    for i in range(len(deltas)):

        val = p_vals[i]

        p_vals[i] -= deltas[i] / 2.
        mp.set_p(p_vals)
        rev_state = mp.get_state()

        p_vals[i] += deltas[i]
        mp.set_p(p_vals)
        fwd_state = mp.get_state()

        fd_grad.append((fwd_state - rev_state) / deltas[i])
        p_vals[i] = val

    # return to the initial state
    mp.set_p(p_vals)

    return fd_grad

def print_model_geometry(beam = None, detector = None, crystal = None):

    if beam:
        print "beam s0 = (%.4f, %.4f, %.4f)" % beam.get_s0()
    if detector:
        print "sensor origin = (%.4f, %.4f, %.4f)" % detector.sensors()[0].origin
        print "sensor dir1 = (%.4f, %.4f, %.4f)" % detector.sensors()[0].dir1
        print "sensor dir2 = (%.4f, %.4f, %.4f)" % detector.sensors()[0].dir2
    if crystal:
        uc = crystal.get_unit_cell()
        print "crystal unit cell = %.4f, %.4f, %.4f, %.4f, %.4f, %.4f" % uc.parameters()
        print "crystal orientation matrix U ="
        print crystal.get_U().round(4)
