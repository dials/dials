# -*- coding: utf-8 -*-

"""This module accepts the data read from diffraction image headers and
computes the coordinates of the crystal and detector using numpy.
"""

import numpy as np
from numpy.linalg import norm

from copy import deepcopy

import viewer

v = dict()        # sample vertices/normals
f = dict()        # sample faces
r_f = dict()      # sample rotated faces
p_v = dict()      # panel vertices/normals
p_f = dict()      # panel faces


class Sample(object):

    """ The crystal class. A cuboid crystal is constructed with the
    origin at its centre, taking the x, y, z dimensions as args.
    """
    def __init__(self, x_length = 4.0, y_length = 0.5, z_length = 8.0):

        self.x_length = x_length
        self.y_length = y_length
        self.z_length = z_length
        x = 0.5 * x_length
        y = 0.5 * y_length
        z = 0.5 * z_length

        # These are initial vectors for each vertex and normal to each
        # face (in each case v[n,0] is the normal vector).
        v[0, 0] = np.array([0.0, 0.0, 1.0])
        v[0, 1] = np.array([x, y, z])
        v[0, 2] = np.array([-x, y, z])
        v[0, 3] = np.array([-x, y, -z])
        v[0, 4] = np.array([x, y, -z])
        # As the cuboid is centred on the origin, opposite faces have
        # -1 * each others' coordinates.
        for i in range(5):
            v[1, i] = -v[0, i]

        v[2, 0] = np.array([0.0, 1.0, 0.0])
        v[2, 1] = np.array([x, y, z])
        v[2, 2] = np.array([x,-y, z])
        v[2, 3] = np.array([-x,-y, z])
        v[2, 4] = np.array([-x, y, z])

        for i in range(5):
            v[3, i] = -v[2, i]

        v[4, 0] = np.array([1.0, 0.0, 0.0])
        v[4, 1] = np.array([x, y, z])
        v[4, 2] = np.array([x, y, -z])
        v[4, 3] = np.array([x,-y, -z])
        v[4, 4] = np.array([x,-y, z])

        for i in range(5):
            v[5, i] = -v[4, i]

        # This loop creates a 2d array for each face, and a dictionary
        # of all faces.
        for p in range(6):
            f[p] = np.zeros((5, 3))
            for q in range(5):
                f[p][q] = v[p,q]
       # Transposing each face array to make it a matrix of column
       # vectors, which allows use of np.dot()
            f[p] = np.transpose(f[p])

        self.f = f

    def rotate(self, i = (1.0, 0.0, 0.0), j = (0.0, 1.0, 0.0),
              k = (0.0, 0.0, 1.0)):
        '''
        Rotates the crystal to a new orientation, given new basis
        vectors as args.
        '''

        sample_basis = np.array([i, j, k])
        # Adjusting each face to its new orientation
        for n in range(6):
        # Transposing again, so that elements are indexed in a way more
        # useful to the viewer.
            r_f[n] = np.transpose(np.dot(sample_basis, self.f[n]))
        self.r_f = r_f
        # r_f is a dictionary containing six 3 by 5 arrays, all the
        # information needed to draw the sample.
        return self.r_f

class Panel(object):

    '''
    The Panel class represents an individual panel of a detector.
    It is constructed using the fast and slow vectors and dimensions,
    and the panel origin.
    '''
    def __init__(self, fast_length = 83.8, slow_length = 106.5, depth = 10.0):
        self.fast_length = fast_length
        self.slow_length = slow_length
        self.depth = depth

    def get_coords(self, origin, fast, slow):

        z = self.depth

        p_v[0, 1] = origin
        p_v[0, 2] = origin + fast * self.fast_length
        p_v[0, 3] = p_v[0, 2] + slow * self.slow_length
        p_v[0, 4] = origin + slow * self.slow_length
        # Cross product of 2 vectors in the plane of the detector face,
        # giving the normal.
        p_v[0, 0] = np.cross(origin - p_v[0, 2], origin - p_v[0, 4])
        p_v[0, 0] /= norm(p_v[0, 0])

        for i in range(1, 5):
            p_v[1, i] = deepcopy(p_v[0, i])
            p_v[1, i] -= z * p_v[0, 0]

        p_v[2, 1] = p_v[0, 2]
        p_v[2, 2] = p_v[0, 3]
        p_v[2, 3] = p_v[1, 3]
        p_v[2, 4] = p_v[1, 2]
        p_v[2, 0] = np.cross(p_v[2, 1] - p_v[2, 2], p_v[2, 1] - p_v[2, 4])
        p_v[2, 0] /= norm(p_v[2, 0])

        for i in range(1, 5):
            p_v[3, i] = deepcopy(p_v[2, i])
            p_v[3, i] += self.fast_length * p_v[2, 0]

        p_v[4, 1] = p_v[0, 1]
        p_v[4, 2] = p_v[1, 1]
        p_v[4, 3] = p_v[1, 2]
        p_v[4, 4] = p_v[0, 2]
        p_v[4, 0] = np.cross(p_v[4, 1] - p_v[4, 2], p_v[4, 1] - p_v[4, 4])
        p_v[4, 0] /= norm(p_v[4, 0])

        for i in range(1,5):
            p_v[5, i] = deepcopy(p_v[4, i])
            p_v[5, i] -= self.slow_length * p_v[4, 0]
        # These loops initialise the normals, create empty face arrays
        # in a dictionary, and then fill them.
        for j in (1, 3, 5):
            p_v[j, 0] = -p_v[j - 1, 0]

        for i in range(6):
            p_f[i] = np.zeros((5, 3))
            for j in range(5):
                p_f[i][j] = p_v[i, j]

        if p_f[0][0][2] != 0.0:
            viewer.DET_ANGLE = -np.rad2deg(np.arctan(p_f[0][0][1]
                                            / p_f[0][0][2]))
        else:
            viewer.DET_ANGLE = np.sign(p_f[0][0][1]) * 90.0
        viewer.detector_normal = p_v[0, 0]
        centre = (origin + 0.5 * fast * self.fast_length
                        + 0.5 * slow * self.slow_length)
        #viewer.d = norm(centre)
        return p_f

if __name__ == "__main__":

    # These are sample i, j, k vectors from a cbf image.
    crystal_i = [ 0.7065882430629097, 0.5534899064860146,-0.44088771607221444]
    crystal_j = [ 0.4514320564090011,-0.8323808627606581,-0.32148281098087234]
    crystal_k = [-0.5449239884714299, 0.02812512627395258,-0.8380136180638496]

    NaCl = Sample()
    r_f = NaCl.rotate(crystal_i, crystal_j, crystal_k)
    print NaCl.r_f, NaCl.r_f[1][2, 0], NaCl.f

    # This is a sample normal from the same data.
    detector_normal = np.array([0.0, -0.4226182617406947, -0.9063077870366522])
    origin = np.array([-41.05, 3.4657241967247856, -116.367397306085])
    fast = np.array([1.0, 0.0, 0.0])
    slow = np.array([0.0, -0.906307787036628, 0.4226182617407468])
    slow_length, fast_length = (83.67757205642768, 106.35843638412338)
    d = Panel()
    d.get_coords(origin, fast, slow)
