from math import exp
import matplotlib.pyplot as plt
import numpy as np
from data_plotter import load_data, plot_data
from data_quality_assessment import R_meas, R_pim
from dials.array_family import flex
from scitbx import lbfgs

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
        print "functional: %12.6g" % f, "gradient norm: %12.6g" % g.norm()
        return f, g

    def residual(self):
        gvalues = self.data_manager.g_values[0:self.data_manager.ndbins]
        resolution = self.res_values
        R = 0.0
        for i, val in enumerate(resolution):
            R += (gvalues[i] * exp((self.x[0])*val) - self.x[1])**2
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


'''
def scale_gvalues(data):
    print list(Optimal_rescale_values.x)
    Optimal_rescale_values = B_optimiser(data, flex.double([0.7, 10.0]))

    scaling_factors = []
    for i in range(0, data.data_manager.nzbins):
        scaling_factors += flex.exp(Optimal_rescale_values.x[0]*Optimal_rescale_values.res_values)
    scaling_factors = flex.double(scaling_factors)

    data.data_manager.g_values = data.data_manager.g_values * scaling_factors
    data.data_manager.g_values = data.data_manager.g_values * (1.0/Optimal_rescale_values.x[1])
    data.data_manager.calc_Ih()
    data.x = data.data_manager.g_values
    print "scaled by Brel and global scale parameter"
    print data.residual()
    return data'''


if __name__ == "__main__":
    filename = "/Users/whi10850/Documents/test_data/integrate/13_integrated_scaled.txt"
    data = load_data(filename)

    scaled_data = scale_gvalues(data)

    Rmeas = R_meas(scaled_data.data_manager)
    Rpim = R_pim(scaled_data.data_manager)
    print "R_meas of the (unmerged) data is %s" % (Rmeas)
    print "R_pim of the merged data is %s" % (Rpim)
    plot_data(scaled_data)
