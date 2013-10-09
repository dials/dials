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
        #tmp_numpy_layering_n_bkgr_modl(reflections)
        layering_and_background_modl(reflections)
        return reflections

def tmp_numpy_layering_n_bkgr_modl(reflections):
    import numpy
    from scitbx.array_family import flex
    print "modelling background tmp numpy"
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

        ref.shoebox = flex.double(shoebox)
        ref.shoebox_background = flex.double(background)

    return reflections

def layering_and_background_modl(reflections):
    from dials.algorithms.background import curved_background_flex_2d
    from scitbx.array_family import flex

    print "performing curved background calculation ...."

    for ref in reflections:
        if ref.is_valid():
            shoebox = ref.shoebox
            mask = ref.shoebox_mask
            background = ref.shoebox_background
            for i in range(shoebox.all()[0]):
                data2d = shoebox[i:i + 1, :, :]
                mask2d = mask[i:i + 1, :, :]
                data2d.reshape(flex.grid(shoebox.all()[1:]))
                mask2d.reshape(flex.grid(shoebox.all()[1:]))
                background2d = curved_background_flex_2d(data2d.as_double(), mask2d)
                background2d.reshape(flex.grid(1, background2d.all()[0], background2d.all()[1]))
                background[i:i + 1, :, :] = background2d.as_double()
    print "curved background calculation .... done"

    return reflections


def layering_and_background_plane(reflections):
    from dials.algorithms.background import plane_background_flex_2d
    from scitbx.array_family import flex

    print "performing background plane calculation ...."

    for ref in reflections:
        if ref.is_valid():
            shoebox = ref.shoebox
            mask = ref.shoebox_mask
            background = ref.shoebox_background
            for i in range(shoebox.all()[0]):
                data2d = shoebox[i:i + 1, :, :]
                mask2d = mask[i:i + 1, :, :]
                data2d.reshape(flex.grid(shoebox.all()[1:]))
                mask2d.reshape(flex.grid(shoebox.all()[1:]))
                background2d = plane_background_flex_2d(data2d.as_double(), mask2d)
                background2d.reshape(flex.grid(1, background2d.all()[0], background2d.all()[1]))
                background[i:i + 1, :, :] = background2d.as_double()
    print "background plane calculation .... done"

    return reflections
