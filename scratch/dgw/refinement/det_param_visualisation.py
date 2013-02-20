#!/usr/bin/env python

# This is a temporary development test script for the detector_parameters
# module. Its purpose is to aid development by providing 3d visualisation
# of geometry relevant to the detector and sensor models, using R and the
# rgl library
from __future__ import division
from detector_parameters import *

def vis(dp, header=True):
    # print out rgl code for testing by visualisation
    if header:
        print "library(rgl)"
        print "colz <- rainbow(10)"
        print "n <- 1"
    print "dorg <- c(%.9f, %.9f, %.9f)" % dp._initial_state['dorg'].elems
    print "offset <- c(%.9f, %.9f, %.9f)" % dp._initial_state['offset'].elems
    print "sensor_origin <- c(%.9f, %.9f, %.9f)" % dp._models[0].origin
    print "sensor_d1 <- c(%.9f, %.9f, %.9f)" % dp._models[0].dir1
    print "sensor_d2 <- c(%.9f, %.9f, %.9f)" % dp._models[0].dir2
    print "sensor_dn <- c(%.9f, %.9f, %.9f)" % dp._models[0].normal
    print "det_d1 <- c(%.9f, %.9f, %.9f)" % dp._initial_state['d1'].elems
    print "det_d2 <- c(%.9f, %.9f, %.9f)" % dp._initial_state['d2'].elems
    print "det_dn <- c(%.9f, %.9f, %.9f)" % dp._initial_state['dn'].elems
    print "sensor_lim1 <- c(%.9f, %.9f)" % dp._models[0].lim1
    print "sensor_lim2 <- c(%.9f, %.9f)" % dp._models[0].lim2
    print 'lines3d(rbind(c(0,0,0),dorg),col="red")'
    print 'lines3d(rbind(c(0,0,0),sensor_origin),col="green")'
    print 'lines3d(rbind(dorg,dorg + offset[1] * det_d1 + offset[2] * det_d2))'
    print "sensor_vertex1 <- sensor_origin + sensor_lim1[1] * sensor_d1 + sensor_lim2[1] * sensor_d2"
    print "sensor_vertex2 <- sensor_origin + sensor_lim1[2] * sensor_d1 + sensor_lim2[1] * sensor_d2"
    print "sensor_vertex3 <- sensor_origin + sensor_lim1[2] * sensor_d1 + sensor_lim2[2] * sensor_d2"
    print "sensor_vertex4 <- sensor_origin + sensor_lim1[1] * sensor_d1 + sensor_lim2[2] * sensor_d2"
    print "sensor_vertices <- matrix(data = c(sensor_vertex1,sensor_vertex2,sensor_vertex3,sensor_vertex4),ncol=3,byrow=T)"
    print "quads3d(sensor_vertices, col=colz[n])"
    print "n <- n + 1"

if __name__ == '__main__':

    # set up a simple detector frame with directions aligned with
    # principal axes and sensor origin located on the z-axis at -100
    d1 = matrix.col((1, 0, 0))
    d2 = matrix.col((0, -1, 0))
    lim = (0,50)
    panel0 = sensor(matrix.col((-20, 20, -100)), d1, d2, lim, lim)
    dp = detector_parameterisation_single_sensor(panel0)

    print "#initial position"
    vis(dp)
    print

    p_vals = dp.get_p()
    p_vals[4:6] = [pi/12, pi/12]
    dp.set_p(dic)
    dp.compose()

    print "#rotation about d1 and d2 by 15 degrees"
    vis(dp, header=False)
    print

    p_vals = [110., 10., 10., pi/12, pi/12, pi/12]
    dp.set_p(p_vals)
    dp.compose()

    print "#shifts by +10 and all rotations by +15 degrees"
    vis(dp, header=False)
    print
