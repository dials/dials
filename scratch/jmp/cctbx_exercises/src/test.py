
from scitbx import matrix

import math

class ReflectionProfileGrid:
    import numpy
    
    def __init__(self, n_reflections, grid_size, step_size, grid_origin=None):

        # If the origin is not set then set it to the centre of the grid
        if grid_origin == None:
            grid_origin = map(lambda x: int(x / 2), grid_size)

        # Calculate the grid shape and allocate memory for the grid array
        grid_shape = (n_reflections,) + grid_size
        self.grid_data = numpy.zeros(shape=grid_shape, dtype=numpy.int32)

        # Set the grid sizes and origin etc
        self.n_reflections = n_reflections
        self.grid_size = grid_size
        self.step_size = step_size
        self.grid_origin = grid_origin



class ReflectionProfile:
    
    def __init__(self, grid, index):
        
        self.grid = grid




def calculate_point_weights(w, wi, ws):
    """Calculate the weight to give each grid point.
    
    If the point w is not within the eight grid points, this function will
    give misleading answers.
    
    :param w: The point coordinate
    :param wi: The 8 surrounding grid coordinates
    :param ws: The step size in e1, e2, e3
    :returns: The weight at each surrounding grid point
    
    """
    return map(lambda wj: abs(1.0 - abs(wj[0] - w[0]) / ws[0]) * 
                          abs(1.0 - abs(wj[1] - w[1]) / ws[1]) * 
                          abs(1.0 - abs(wj[2] - w[2]) / ws[2]), wi)



def test_calculate_point_weights():

    from random import uniform
    import math
    import numpy

    # Pick random step sizes between 0.01 and 10
    ws = matrix.col((uniform(0.01, 10), 
                     uniform(0.01, 10), 
                     uniform(0.01, 10)))

    # start as zero and generate all the other points
    wi = []
    for k in (0, ws[2]):
        for j in (0, ws[1]):
            for i in (0, ws[0]):
                wi.append(matrix.col((i, j, k)))
    wi = numpy.array(wi)
    
    # pick 100 random points within the points
    wp = []
    for i in range(0, 100):
        wp.append(matrix.col((uniform(0.0, ws[0]),
                              uniform(0.0, ws[1]),
                              uniform(0.0, ws[2]))))

    # Loop through all the random points and calculate the probability
    pi = []
    for w in wp:
        pi.append(calculate_point_weights(w, wi, ws))

    # Check total probability equal to 1
    epsilon = 1e-10
    for p in pi:
        if 1.0 + epsilon < math.fsum(p) < 1.0 - epsilon:
            raise "Sum of pi != 1"
        
    # Check centroid is at w
    for w, p, in zip(wp, pi):
        centroid = matrix.col((math.fsum(p * wi[:,0]), 
                               math.fsum(p * wi[:,1]), 
                               math.fsum(p * wi[:,2])))
        diff = (w - centroid).length()
        if 1.0 + epsilon < diff < 1.0 - epsilon:
            raise "Centroid not at w"

    print "Succeeded"


test_calculate_point_weights()

#def prob(l, wi):
#    
#    elwj = map(lambda wj: math.exp(-l.dot(wj)), wi)
#    sum_elwj = math.fsum(elwj)
#    return map(lambda wj: math.exp(-l.dot(wj)) / sum_elwj, wi)
#
#def prob2(w, wi, ws):
#
##    return map(lambda wj: (abs(ws[0] - abs(wj[0] - w[0])) / ws[0]) * 
##                          (abs(ws[1] - abs(wj[1] - w[1])) / ws[1]) * 
##                          (abs(ws[2] - abs(wj[2] - w[2])) / ws[2]), wi)
#    return map(lambda wj: abs(1.0 - abs(wj[0] - w[0]) / ws[0]) * 
#                          abs(1.0 - abs(wj[1] - w[1]) / ws[1]) * 
#                          abs(1.0 - abs(wj[2] - w[2]) / ws[2]), wi)
#
#ws = [1, 2, 3]
#wi = []
#for k in (0, 3):
#    for j in (0, 2):
#        for i in (0, 1):
#            wi.append(matrix.col((i, j, k)))
            
            
##w = matrix.col((0.1, 0.1, 0.1))
#
#wc = matrix.col((0.5, 1.0, 1.5))
#w = matrix.col((0.67, 0.2, 1.5))
##w = wc
#l = (wc - w)
#l = matrix.col((3, 2, 1))-w
#
##l = w
#
#pi = prob2(w, wi, ws)
#
##pi = prob(l, wi)
##
#for w, p in zip(wi, pi):
#    print tuple(w), p
# #   
#import numpy
#pi = numpy.array(pi)
#wi = numpy.array(wi)
# 
#print pi
#print math.fsum(pi)
#print math.fsum(pi * wi[:,0]), math.fsum(pi * wi[:,1]), math.fsum(pi * wi[:,2])
#


#w0 = matrix.col((0, 0, 0))
#w1 = matrix.col((1, 0, 0))
#w = matrix.col((0.5, 0, 0))
#
#l = w0.cross(w1)
#print l.dot(w0), l.dot(w1)