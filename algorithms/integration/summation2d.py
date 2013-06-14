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
    def __call__(self, sweep, crystal, reflections):
        '''Process the reflections.'''
        self.integrate(reflections)
        return reflections

    def integrate(self, reflections):
        tmp_translating_n_layering(reflections)

def tmp_translating_n_layering(reflections):
    print "hi there translating and layering 01"

    from dials.algorithms.integration.sumation_2d import raw_2d_integration
    for ref in reflections:
        if ref.status == 0:
            shoebox = ref.shoebox.as_numpy_array()
            mask = ref.shoebox_mask.as_numpy_array()
            backgound = ref.shoebox_background.as_numpy_array()

            for i in range(shoebox.shape[0]):
                data2d = shoebox[i]
                mask2d = mask[i]
                backgound2d = backgound[i]
                itns, sigma = raw_2d_integration(data2d, mask2d, backgound2d)
            ref.intensity = float(itns)
            ref.intensity_variance = float(sigma * sigma)

    print "hi there translating and layering 02"
