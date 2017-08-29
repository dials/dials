from dials_array_family_flex_ext import *
from cctbx.array_family.flex import *
from cctbx.array_family import flex
from scitbx import lbfgs
from math import exp


class optimiser(object):
    '''Class that takes in Data_Manager object and runs
    an LBFGS minimisation in a Kabsch approach '''
    def __init__(self, Data_Manager_object, sigma):
        self.data_manager = Data_Manager_object
        self.x = self.data_manager.g_values
        self.sigma = sigma
        self.residuals = []
        lbfgs.run(target_evaluator=self)

    def compute_functional_and_gradients(self):
        self.data_manager.calc_Ih()
        f = self.residual()
        g = self.gradient()
        self.residuals.append(f)
        print "functional: %12.6g" % f, "gradient norm: %12.6g" % g.norm()
        return f, g

    def residual(self):
        intensities = self.data_manager.sorted_reflections['intensity.sum.value']
        g_values = self.x
        variances = self.data_manager.sorted_reflections['intensity.sum.variance']
        Ih_array = self.data_manager.Ih_array
        R = 0.0
        for n in range(len(intensities)):
            l = self.data_manager.sorted_reflections['l_bin_index'][n]
            h = self.data_manager.sorted_reflections['h_index'][n]
            R += (((intensities[n]-(g_values[l]*Ih_array[h]))**2)/variances[n])
        for l in range(len(g_values)):
            R += (((g_values[l]-1)**2)/(self.sigma**2))
        return R

    def gradient(self):
        h_index_counter_array = self.data_manager.h_index_counter_array
        h_index_cumulative_array = self.data_manager.h_index_cumulative_array
        intensities = self.data_manager.sorted_reflections['intensity.sum.value']
        g_values = self.x
        variances = self.data_manager.sorted_reflections['intensity.sum.variance']
        Ih_array = self.data_manager.Ih_array
        G = flex.double([0.0]*len(g_values))
        for n in range(len(h_index_counter_array)):
            lsum = h_index_counter_array[n]
            sumgsq = 0.0
            for i in range(lsum):
                indexer = i + h_index_cumulative_array[n]
                l = self.data_manager.sorted_reflections['l_bin_index'][indexer]
                sumgsq += ((g_values[l]**2)/variances[indexer])
            for i in range(lsum):
                '''determine index of data table'''
                indexer = i + h_index_cumulative_array[n]
                l = self.data_manager.sorted_reflections['l_bin_index'][indexer]
                h = self.data_manager.sorted_reflections['h_index'][indexer]
                dIhdg = ((intensities[indexer]/variances[indexer])
                         - (Ih_array[h]*2.0*g_values[l]/variances[indexer]))/sumgsq
                rhl = (intensities[indexer]-(g_values[l]*Ih_array[h]))/(variances[indexer]**0.5)
                G[l] += 2.0*rhl*((-1.0*Ih_array[h]/(variances[indexer]**0.5))
                                 -((g_values[l]/(variances[indexer]**0.5)) * dIhdg)
                                )
        for l in range(len(G)):
            G[l] += (2.0*(g_values[l]-1.0))/(self.sigma**2.0)
        return G









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
            R += ((gvalues[i] * exp((self.x[0])*val)) - self.x[1])**2
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