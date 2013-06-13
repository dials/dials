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
class Integrate2d(IntegrationInterface):
    '''A class to perform 2D integration'''

    def __init__(self, **kwargs):
        '''Initialise algorithm.'''
        pass
    def __call__(self, reflections):
        '''Process the reflections.'''
        self.subract_background(reflections)
        self.integrate(reflections)
    def subract_background(self, reflections):
        from dials.algorithms.background.background_subtraction_2d \
          import flat_background_calc_2d, curved_background_calc_2d

        import numpy


        from scitbx.array_family import flex
        for ref in reflections:

            shoebox = ref.shoebox.as_numpy_array()
            mask = ref.shoebox_mask.as_numpy_array()
            background = numpy.copy(shoebox)
            for i in range(shoebox.shape[0]):
                data2d = shoebox[i]
                mask2d = mask[i]
                background2d = background[i]

                background2d = flat_background_calc_2d(data2d, mask2d)
                #background2d = curved_background_calc_2d(data2d, mask2d)
                background[i] = background2d
            ref.shoebox = flex.int(shoebox)
            ref.shoebox_background = flex.int(background)

    def integrate(self, reflections):
        from dials.algorithms.integration.sumation_2d import raw_2d_integration
        for ref in reflections:
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
