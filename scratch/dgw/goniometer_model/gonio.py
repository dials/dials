from __future__ import division
from scitbx import matrix

class goniometer:
    '''Dumb goniometer: merely a unit rotation axis'''

    def __init__(self, axis):

        assert isinstance(axis, matrix.rec)
        self._axis = axis.normalize()

    def get_axis(self):
        return self._axis

    def set_axis(self, vals):
        if len(vals) != 3:
            assert isinstance(vals, matrix.col)
        self._axis = matrix.col(vals).normalize()
