import numpy as np
from dials_array_family_flex_ext import *
from cctbx.array_family.flex import *
from cctbx.array_family import flex
from scitbx import lbfgs
from math import exp
#note: include math exp import after flex imports to avoid exp conflicts

class optimiser(object):
    '''Class that takes in Data_Manager object and runs
    an LBFGS minimisation in a Kabsch approach '''
    def __init__(self, Data_Manager_object, correction,
                 decay_correction_rescaling=False):
        self.data_manager = Data_Manager_object
        self.residuals = []
        self.correction = correction
        '''This next part could be condensed using class inheritance'''
        if correction == 'decay':
            print "performing decay correction"
            self.bin_index = 'l_bin_index'
            self.size = (self.data_manager.binning_parameters['nzbins'] * 
                         self.data_manager.binning_parameters['ndbins'])
            self.x = self.data_manager.g_values
            self.gvalues_stat1 = flex.double([self.data_manager.g2_values[i]
                for i in self.data_manager.sorted_reflections['a_bin_index']])
            self.gvalues_stat2 = flex.double([self.data_manager.g3_values[i]
                for i in self.data_manager.sorted_reflections['xy_bin_index']])  
        elif correction == 'absorption':
            print "performing absorption correction"
            self.bin_index = 'a_bin_index'
            self.size = (self.data_manager.binning_parameters['n_absorption_positions'])**2
            self.x = self.data_manager.g2_values
            self.gvalues_stat1 = flex.double([self.data_manager.g_values[i]
                for i in self.data_manager.sorted_reflections['l_bin_index']])  
            self.gvalues_stat2 = flex.double([self.data_manager.g3_values[i]
                for i in self.data_manager.sorted_reflections['xy_bin_index']])
        elif correction == 'modulation':
            print "performing modulation correction"
            self.bin_index = 'xy_bin_index'
            self.size = (self.data_manager.binning_parameters['n_detector_bins'])**2
            self.x = self.data_manager.g3_values
            self.gvalues_stat1 = flex.double([self.data_manager.g_values[i]
                for i in self.data_manager.sorted_reflections['l_bin_index']])
            self.gvalues_stat2 = flex.double([self.data_manager.g2_values[i]
                for i in self.data_manager.sorted_reflections['a_bin_index']])
        else:
            raise ValueError('Invalid bin index type')
        self.core_parameters = lbfgs.core_parameters(maxfev=10)
        self.termination_params = lbfgs.termination_parameters(max_iterations=10)
        lbfgs.run(target_evaluator=self, termination_params=self.termination_params,
                  core_params=self.core_parameters)
        if correction == 'decay':
            if decay_correction_rescaling:
                self.data_manager.scale_gvalues()


    def compute_functional_and_gradients(self):
        '''first calculate the updated values of the necessary arrays'''
        gxvalues = flex.double([self.x[i] for i in
                   self.data_manager.sorted_reflections[self.bin_index]])
        self.data_manager.scale_factors = gxvalues * self.gvalues_stat1 * self.gvalues_stat2
        self.data_manager.calc_Ih()
        f = self.residual()
        g = self.gradient()
        self.residuals.append(f)
        print "Residual: %12.6g" % f
        return f, g

    def residual(self):
        intensities = self.data_manager.sorted_reflections[self.data_manager.int_method[0]]
        variances = self.data_manager.sorted_reflections[self.data_manager.int_method[1]]
        R = (((intensities - (self.data_manager.scale_factors * self.data_manager.Ih_values))**2
             ) / variances)
        R = flex.sum(R)#/1000.0
        return R

    def gradient(self):
        intensities = self.data_manager.sorted_reflections[self.data_manager.int_method[0]]
        variances = self.data_manager.sorted_reflections[self.data_manager.int_method[1]]
        gsq = ((self.data_manager.scale_factors)**2) / variances
        sumgsq = np.add.reduceat(gsq, self.data_manager.h_index_cumulative_array[:-1])
        sumgsq = flex.double(np.repeat(sumgsq, self.data_manager.h_index_counter_array))
        dIhdg = flex.double(((intensities / variances)
                             - (self.data_manager.Ih_values * 2.0
                                * self.data_manager.scale_factors / variances)
                            ) / sumgsq)
        rhl = flex.double((intensities - (self.data_manager.Ih_values
                                          * self.data_manager.scale_factors)
                          ) / (variances**0.5))
        grad = 2.0 * rhl * ((-1.0 * self.data_manager.Ih_values / (variances**0.5))
                            - ((self.data_manager.scale_factors / (variances**0.5)) * dIhdg))
        return flex.double(np.bincount(self.data_manager.sorted_reflections[self.bin_index],
                                       weights=grad, minlength=self.size))

class B_optimiser(object):
    '''Class that takes in Data_Manager object and runs
    an LBFGS minimisation in a Kabsch approach '''
    def __init__(self, Data_Manager_object, initial_values):
        self.data_manager = Data_Manager_object
        d_bin_boundaries = self.data_manager.bin_boundaries['d']
        self.res_values = flex.double([])
        for i in range(0, len(d_bin_boundaries) - 1):
            self.res_values.append(((1.0 / (d_bin_boundaries[i]**2))
                                    +(1.0 / (d_bin_boundaries[i+1]**2))) / 2.0)
        self.x = initial_values
        lbfgs.run(target_evaluator=self)

    def compute_functional_and_gradients(self):
        f = self.residual()
        g = self.gradient()
        return f, g

    def residual(self):
        gvalues = self.data_manager.g_values[0:self.data_manager.binning_parameters['ndbins']]
        resolution = self.res_values
        R = 0.0
        for i, val in enumerate(resolution):
            R += ((gvalues[i] * exp(self.x[0] * val)) - self.x[1])**2
        return R

    def gradient(self):
        gvalues = self.data_manager.g_values[0:self.data_manager.binning_parameters['ndbins']]
        resolution = self.res_values
        G = flex.double([0.0, 0.0])
        for i, val in enumerate(resolution):
            G[0] += (2.0 * ((gvalues[i] * exp((self.x[0]) * val)) - self.x[1])
                     * resolution[i]*gvalues[i]*exp((self.x[0])*resolution[i]))
            G[1] += -2.0 * ((gvalues[i] * exp((self.x[0]) * val)) - self.x[1])
        return G
