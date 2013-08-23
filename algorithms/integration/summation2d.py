
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
    def __call__(self, sweep, crystal, reflections, reference=None):
        '''Process the reflections.'''
        self.integrate(reflections)
        return reflections

    def integrate(self, reflections):
        #tmp_numpy_layering_n_integrating(reflections)
        flex_2d_layering_n_integrating(reflections)

def tmp_numpy_layering_n_integrating(reflections):
    print "hi there translating and layering 01 "

    from dials.algorithms.integration.summation_2d import raw_2d_integration
    for ref in reflections:
        if ref.is_valid():
            shoebox = ref.shoebox.as_numpy_array()
            mask = ref.shoebox_mask.as_numpy_array()
            backgound = ref.shoebox_background.as_numpy_array()
            ref.intensity = 0.0
            ref.intensity_variance = 0.0
            for i in range(shoebox.shape[0]):
                data2d = shoebox[i]
                mask2d = mask[i]
                backgound2d = backgound[i]
                itns, sigma = raw_2d_integration(data2d, mask2d, backgound2d)
                ref.intensity += float(itns)
                ref.intensity_variance += float(sigma * sigma)
                #ref.intensity_variance += float(sigma)
            #ref.intensity_variance = ref.intensity_variance * ref.intensity_variance

    print "hi there translating and layering 02 "
def flex_2d_layering_n_integrating(reflections):
    from scitbx.array_family import flex
    from dials.algorithms.integration import raw_2d_cut
    print "starts flex array version of summation integration "

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
                ref.intensity_variance += reslt[1] * reslt[1]

    print "flex array version of summation integration ends"

    return reflections
