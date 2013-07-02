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
        tmp_numpy_layering_n_integrating(reflections)

def tmp_numpy_layering_n_integrating(reflections):
    print "hi there translating and layering 01"

    from dials.algorithms.integration.summation_2d import raw_2d_integration
    for ref in reflections:
        if ref.is_valid():
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
def flex_2d_layering_n_integrating(reflections):
    from scitbx.array_family import flex
    from dials.scratch.luiso_s import raw_2d_cut
    print "flex array version of summation integration 01"

    for ref in reflections:
        if ref.is_valid():
            shoebox = ref.shoebox
            mask = ref.shoebox_mask
            background = ref.shoebox_background

            for i in range(shoebox.all()[0]):
                data2d = shoebox[i:i + 1, :, :]
                mask2d = mask[i:i + 1, :, :]
                background2d = background[i:i + 1, :, :]

                data2d.reshape(flex.grid(shoebox.all()[1:]))
                mask2d.reshape(flex.grid(shoebox.all()[1:]))
                background2d.reshape(flex.grid(shoebox.all()[1:]))

                reslt = raw_2d_cut(data2d, mask2d, background2d)

            ref.intensity = reslt[0]
            ref.intensity_variance = reslt[1]

    return reflections

##############################################################################################
    from dials.algorithms.integration.summation_2d import raw_2d_integration
    for ref in reflections:
        if ref.is_valid():
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
#############################################################################################
    a = tst_01(flex.int(data2d)).as_numpy_array()
    print a
    '''
        from dials.algorithms.background import curved_background_flex_2d
        from scitbx.array_family import flex
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
    
    
        return reflections
    '''
    print "flex array version of summation integration 02"
