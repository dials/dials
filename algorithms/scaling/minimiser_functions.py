'''
Classes to create minimiser objects.
'''

import numpy as np
from dials_array_family_flex_ext import *
from cctbx.array_family.flex import *
from cctbx.array_family import flex
from scitbx import lbfgs
from math import exp
from target_function import *
#note: include math exp import after flex imports to avoid exp conflicts

class LBFGS_optimiser(object):
    '''Class that takes in Data_Manager object and runs
    an LBFGS minimisation in a Kabsch approach '''
    def __init__(self, Data_Manager_object, parameters, param_name, 
                 decay_correction_rescaling=False):
        self.data_manager = Data_Manager_object
        self.x = parameters
        self.parameter_name = param_name
        self.set_up_parameterisation()
        self.residuals = []
        lbfgs.run(target_evaluator=self)
        if decay_correction_rescaling:
            self.data_manager.scale_gvalues()

    def set_up_parameterisation(self):
        '''Set up the problem by indicating which g values are being minimised'''
        constant_g_values = []
        for parameterisation_type in self.data_manager.sf_parameterisations:
            index = self.data_manager.sf_parameterisations.index(parameterisation_type)
            bin_index = self.data_manager.bin_indices[index]
            if self.parameter_name == parameterisation_type:
                #change name to active bin index?
                self.data_manager.bin_index = bin_index
                self.data_manager.param_size = len(self.x)
            else:
                constant_g_values.append(
                    flex.double([self.data_manager.g_parameterisation[index][i]
                                 for i in self.data_manager.sorted_reflections[bin_index]]))
        self.data_manager.constant_g_values = constant_g_values[0]*constant_g_values[1]

    def compute_functional_and_gradients(self):
        '''first calculate the updated values of the scale factors and Ih,
        before calculating the residual and gradient functions'''
        self.data_manager.update_for_minimisation(parameters=self.x)
        f, g = self.data_manager.get_target_function()
        f = flex.sum(f)
        self.residuals.append(f)
        print "Residual sum: %12.6g" % f
        return f, g


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
        gvalues = self.data_manager.g_decay[0:self.data_manager.binning_parameters['ndbins']]
        resolution = self.res_values
        R = 0.0
        for i, val in enumerate(resolution):
            R += ((gvalues[i] * exp(self.x[0] * val)) - self.x[1])**2
        return R

    def gradient(self):
        gvalues = self.data_manager.g_decay[0:self.data_manager.binning_parameters['ndbins']]
        resolution = self.res_values
        G = flex.double([0.0, 0.0])
        for i, val in enumerate(resolution):
            G[0] += (2.0 * ((gvalues[i] * exp((self.x[0]) * val)) - self.x[1])
                     * resolution[i]*gvalues[i]*exp((self.x[0])*resolution[i]))
            G[1] += -2.0 * ((gvalues[i] * exp((self.x[0]) * val)) - self.x[1])
        return G
