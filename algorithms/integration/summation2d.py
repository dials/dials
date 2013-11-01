
#
# integrate2d.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: Luiso & James
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dials.interfaces.integration import IntegrationInterface
class Summation2d(IntegrationInterface):
    '''A class to perform 2D integration'''

    def __init__(self, **kwargs):
        '''Initialise algorithm.'''
        pass
    def __call__(self, sweep, crystal, reflections, reference = None):
        '''Process the reflections.'''
        self.integrate(reflections)
        return reflections

    def integrate(self, reflections):
        flex_2d_layering_n_integrating(reflections)

def flex_2d_layering_n_integrating(reflections):
    from scitbx.array_family import flex
    from dials.algorithms.integration import raw_2d_cut
    print "performing summation integration .... "

    for ref in reflections:
        if ref.is_valid():
            shoebox = ref.shoebox
            mask = ref.shoebox_mask
            background = ref.shoebox_background
            ref.intensity = 0.0
            ref.intensity_variance = 0.0
            for i in range(shoebox.all()[0]):
                data2d = shoebox[i:i + 1, :, :]
                mask2d = mask[i:i + 1, :, :]
                background2d = background[i:i + 1, :, :]

                data2d.reshape(flex.grid(shoebox.all()[1:]))
                mask2d.reshape(flex.grid(shoebox.all()[1:]))
                background2d.reshape(flex.grid(shoebox.all()[1:]))

                reslt = raw_2d_cut(data2d, mask2d, background2d)

                ref.intensity += reslt[0]
                ref.intensity_variance += reslt[1]

    print "summation integration .... done"

    return reflections
