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
n_ref = 0
ref_bkgr = []
class Integrate2d(IntegrationInterface):
    '''A class to perform 2D integration'''

    def __init__(self, **kwargs):
        '''Initialise algorithm.'''
        pass
    def __call__(self, reflections):
        '''Process the reflections.'''
        print 'n_ref =', n_ref
        print 'len(ref_bkgr) =', len(ref_bkgr)
        self.subract_background(reflections)
        print 'n_ref =', n_ref
        print 'len(ref_bkgr) =', len(ref_bkgr)
        self.integrate(reflections)
    def subract_background(self, reflections):
        global n_ref, ref_bkgr
        from dials.algorithms.background.background_subtraction_2d \
          import flat_background_subtraction_2d , curved_background_subtraction_2d, flat_background_calc_2d
        from dials.algorithms.background import background_subtract_2d

        import numpy


        from scitbx.array_family import flex
        n_ref = len(reflections)
        for ref in reflections:

            shoebox = ref.shoebox.as_numpy_array()
            mask = ref.shoebox_mask.as_numpy_array()
            background = numpy.copy(shoebox)
            #tot_bkgr = 0
            #cont_bkgr = 0
            for i in range(shoebox.shape[0]):
                data2d = shoebox[i]
                mask2d = mask[i]
                background2d = background[i]
                ######################################################################################
                #fl_Bkg = background_subtract_2d(flex.int(data2d))
                #data2d = fl_Bkg.as_numpy_array()
                #bkgr = 2.625
                ######################################################################################

                #bkgr = flat_background_calc_2d(data2d, mask2d)

                background2d = flat_background_calc_2d(data2d, mask2d)


                #bkgr = flat_background_subtraction_2d(data2d, mask2d)
                #bkgr = curved_background_subtraction_2d(data2d, mask2d)

                #tot_bkgr += 2.625
                #cont_bkgr += 1
                #shoebox[i] = data2d
                background[i] = background2d
            #avg_bkgr = tot_bkgr / cont_bkgr
            #ref_bkgr.append(avg_bkgr)
            #print 'avg_bkgr =', avg_bkgr
            ref.shoebox = flex.int(shoebox)
            ref.shoebox_background = flex.int(background)

    def integrate(self, reflections):
        from dials.algorithms.integration.sumation_2d import raw_2d_integration
        global ref_bkgr
        cont = 0
        for ref in reflections:
            shoebox = ref.shoebox.as_numpy_array()
            mask = ref.shoebox_mask.as_numpy_array()
            backgound = ref.shoebox_background.as_numpy_array()

            for i in range(shoebox.shape[0]):
                data2d = shoebox[i]
                mask2d = mask[i]
                backgound2d = backgound[i]
                #itns, sigma = raw_2d_integration(data2d, mask2d, ref_bkgr[cont])
                itns, sigma = raw_2d_integration(data2d, mask2d, backgound2d)
            cont += 1
            ref.intensity = float(itns)
            ref.intensity_variance = float(sigma * sigma)
