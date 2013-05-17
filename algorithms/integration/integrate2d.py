#
# integrate2d.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
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

#        from dials.algorithms.background.flat_subtraction import flat_background_subtraction_2d # <- here
#        from dials.algorithms.background.curved_subtraction import curved_background_subtraction_2d # <- here
        #from dials.algorithms.background.background_subtraction_2d  import curved_background_subtraction_2d # <- here

        from dials.algorithms.background.background_subtraction_2d \
          import flat_background_subtraction_2d , curved_background_subtraction_2d

        from scitbx.array_family import flex
        for ref in reflections:
            shoebox = ref.shoebox.as_numpy_array()
            mask = ref.shoebox_mask.as_numpy_array()

            for i in range(shoebox.shape[0]):
                data2d = shoebox[i]
                mask2d = mask[i]

                bkgr = flat_background_subtraction_2d(data2d, mask2d) # <- and here
#                bkgr = curved_background_subtraction_2d(data2d, mask2d) # <- and here

                shoebox[i] = data2d
                mask[i] = mask2d

            ref.shoebox = flex.int(shoebox)
            ref.mask = flex.int(mask)

    def integrate(self, reflections):
        from dials.algorithms.integration.sumation_2d import raw_2d_integration

        for ref in reflections:
            shoebox = ref.shoebox.as_numpy_array()
            mask = ref.shoebox_mask.as_numpy_array()

            for i in range(shoebox.shape[0]):
                data2d = shoebox[i]
                mask2d = mask[i]
                itns, sigma = raw_2d_integration(data2d, mask2d)
            ref.intensity = float(itns)
            ref.intensity_variance = float(sigma)
