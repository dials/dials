#!/usr/bin/env python

"""This model is deprecated. The code remains here only to show how ideas
developed over time"""

# fixme boilerplate
from __future__ import division

from scitbx import matrix

class sensor:
    '''A little class to represent one sensor, which is assumed to be a
    flat rectangle defined by: an origin, two directions and extents in those
    two directions. This may or may not map onto the pixel directions on the
    image. All parameters are cctbx vectors or tuples.'''

    def __init__(self, origin, dir1, dir2, lim1, lim2):
        self._origin = origin
        self._dir1 = dir1
        self._dir2 = dir2
        self._lim1 = lim1
        self._lim2 = lim2

        self._normal = dir1.cross(dir2)
        self._distance = origin.dot(self._normal)

        return

    def intersect(self, ray):
        '''Compute intersection of sensor with ray from frame origin, returning
        none if intersection not within limits.'''
        scale = ray.dot(self._normal)
        r = (ray * self._distance / scale) - self._origin
        x1 = r.dot(self._dir1)
        x2 = r.dot(self._dir2)
        if x1 < self._lim1[0] or x1 > self._lim1[1]:
            return None
        if x2 < self._lim2[0] or x2 > self._lim2[1]:
            return None
        return (x1, x2)

def test_sensor():
    origin = matrix.col((-10, -10, 100))
    dir1 = matrix.col((1, 0, 0))
    dir2 = matrix.col((0, 1, 0))
    lim1 = (0, 20)
    lim2 = (0, 20)

    s = sensor(origin, dir1, dir2, lim1, lim2)

    r = matrix.col((0, 0, 1))
    assert(s.intersect(r))

    r = matrix.col((0, 1, 1))
    assert(not s.intersect(r))

    return

class detector:
    '''An abstract class definition for X-ray detectors which are assumed to
    be composed of one or more flat rectangular sensor areas. Will initially
    assume that sensor areas to not shadow one another.'''

    def __init__(self, sensors):
        self._sensors = sensors

        return

if __name__ == '__main__':
    test_sensor()
