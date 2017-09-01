import numpy as np
from dials_array_family_flex_ext import *
from cctbx.array_family.flex import *
from cctbx.array_family import flex
from scitbx import lbfgs
from math import exp

class optimiser(object):
    '''Class that takes in Data_Manager object and runs
    an LBFGS minimisation in a Kabsch approach '''
    def __init__(self, Data_Manager_object, sigma, correction):
        self.data_manager = Data_Manager_object
        self.sigma = sigma
        self.residuals = []
        self.correction = correction
        '''This next part could be condensed using class inheritance'''
        if correction == 'decay':
            print "performing decay correction"
            self.bin_index = 'l_bin_index'
            self.size = (self.data_manager.nzbins * self.data_manager.ndbins)
            self.x = self.data_manager.g_values
            self.gvalues_stat1 = flex.double([self.data_manager.g2_values[i] for i in self.data_manager.sorted_reflections['a_bin_index']])
            self.gvalues_stat2 = flex.double([self.data_manager.g3_values[i] for i in self.data_manager.sorted_reflections['xy_bin_index']])  
        elif correction == 'absorption':
            print "performing absorption correction"
            self.bin_index = 'a_bin_index'
            self.size = (self.data_manager.npos)**2
            self.x = self.data_manager.g2_values
            self.gvalues_stat1 = flex.double([self.data_manager.g_values[i] for i in self.data_manager.sorted_reflections['l_bin_index']])  
            self.gvalues_stat2 = flex.double([self.data_manager.g3_values[i] for i in self.data_manager.sorted_reflections['xy_bin_index']])
        elif correction == 'modulation':
            print "performing modulation correction"
            self.bin_index = 'xy_bin_index'
            self.size = (self.data_manager.ngridpoints)**2
            self.x = self.data_manager.g3_values
            self.gvalues_stat1 = flex.double([self.data_manager.g_values[i] for i in self.data_manager.sorted_reflections['l_bin_index']])
            self.gvalues_stat2 = flex.double([self.data_manager.g2_values[i] for i in self.data_manager.sorted_reflections['a_bin_index']])
        else:
            raise ValueError('Invalid bin index type')
        self.core_parameters = lbfgs.core_parameters(maxfev=10)
        self.termination_params = lbfgs.termination_parameters(max_iterations=10)
        lbfgs.run(target_evaluator=self, termination_params=self.termination_params, core_params=self.core_parameters)

    def compute_functional_and_gradients(self):
        '''first calculate the updated values of the necessary arrays'''
        self.data_manager.calc_Ih()
        gxvalues = flex.double([self.x[i] for i in self.data_manager.sorted_reflections[self.bin_index]]) 
        self.Ihvalues = flex.double([self.data_manager.Ih_array[i] for i in self.data_manager.sorted_reflections['h_index']])
        self.gproduct = gxvalues * self.gvalues_stat1 * self.gvalues_stat2
        f = self.residual()
        g = self.gradient()
        self.residuals.append(f)
        print "Residual: %12.6g" % f
        return f, g

    def residual(self):    
        intensities = self.data_manager.sorted_reflections['intensity.sum.value']
        variances = self.data_manager.sorted_reflections['intensity.sum.variance']
        R = (((intensities - (self.gproduct * self.Ihvalues))**2)/variances)
        R = np.sum(R)/1000.0
        return R

    def gradient(self):
        h_index_counter_array = self.data_manager.h_index_counter_array
        intensities = self.data_manager.sorted_reflections['intensity.sum.value']
        variances = self.data_manager.sorted_reflections['intensity.sum.variance']

        gsq = ((self.gproduct)**2)/variances
        sumgsq = np.bincount(self.data_manager.sorted_reflections['h_index'], gsq)
        sumgsq = np.repeat(sumgsq, h_index_counter_array)

        dIhdg = ((intensities/variances) - (self.Ihvalues * 2.0 * self.gproduct/variances))/sumgsq
        rhl = (intensities - (self.Ihvalues * self.gproduct))/(variances**0.5)
        grad = 2.0*rhl*((-1.0*self.Ihvalues/(variances**0.5))- ((self.gproduct/(variances**0.5))*dIhdg))
    
        return flex.double(np.bincount(self.data_manager.sorted_reflections[self.bin_index], weights=grad, minlength=self.size))
        



class B_optimiser(object):
    '''Class that takes in Data_Manager object and runs
    an LBFGS minimisation in a Kabsch approach '''
    def __init__(self, Data_Manager_object, initial_values):
        self.data_manager = Data_Manager_object
        d_bin_boundaries = self.data_manager.bin_boundaries['d']
        self.res_values = flex.double([])
        for i in range(0, len(d_bin_boundaries)-1):
            self.res_values.append(((1.0/(d_bin_boundaries[i]**2))
                                    +(1.0/(d_bin_boundaries[i+1]**2)))/2.0)
        self.x = initial_values
        lbfgs.run(target_evaluator=self)

    def compute_functional_and_gradients(self):
        f = self.residual()
        g = self.gradient()
        return f, g

    def residual(self):
        gvalues = self.data_manager.g_values[0:self.data_manager.ndbins]
        resolution = self.res_values
        R = 0.0
        for i, val in enumerate(resolution):
            R += ((gvalues[i] * exp(self.x[0]*val)) - self.x[1])**2
        return R

    def gradient(self):
        gvalues = self.data_manager.g_values[0:self.data_manager.ndbins]
        resolution = self.res_values
        G = flex.double([0.0, 0.0])
        for i, val in enumerate(resolution):
            G[0] += (2.0 * ((gvalues[i] * exp((self.x[0])*val)) - self.x[1])
                     * resolution[i]*gvalues[i]*exp((self.x[0])*resolution[i]))
            G[1] += -2.0 * ((gvalues[i] * exp((self.x[0])*val)) - self.x[1])
        return G