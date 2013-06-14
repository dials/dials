#
# dials.algorithms.background.curved_subtractor.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: Luis Fuentes-Montero "luiso" & James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dials.interfaces.background import BackgroundSubtractionInterface

from dials.algorithms.background.background_subtraction_2d \
          import curved_background_calc_2d

class CurvedSubtractor(BackgroundSubtractionInterface):
    ''' The Flat background subtractor '''

    def __init__(self, **kwargs):
        pass

    def __call__(self, sweep, crystal, reflections):
        import numpy

        from scitbx.array_family import flex
        print "modeling background"

        for ref in reflections:
            shoebox = ref.shoebox.as_numpy_array()
            mask = ref.shoebox_mask.as_numpy_array()
            background = numpy.copy(shoebox)
            for i in range(shoebox.shape[0]):
                data2d = shoebox[i]
                mask2d = mask[i]
                background2d = background[i]

                background2d = curved_background_calc_2d(data2d, mask2d)
                background[i] = background2d
            ref.shoebox = flex.int(shoebox)
            ref.shoebox_background = flex.int(background)
        # Return the reflections
        return reflections
