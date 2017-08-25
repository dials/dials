from dials_array_family_flex_ext import *
from cctbx.array_family.flex import *
from cctbx.array_family import flex
from scitbx import lbfgs


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











def calc_dIhdg(data_manager, Ih_array):
    n_bins = data_manager.n_bins
    n_unique_indices = data_manager.n_unique_indices
    intensities = data_manager.sorted_reflections['intensity.sum.value']
    g_values = data_manager.g_values
    variances = data_manager.sorted_reflections['intensity.sum.variance']

    dIh_array = flex.double([0.0] * n_bins * n_unique_indices)
    dIh_array.reshape(flex.grid(n_bins, n_unique_indices))
    numerator = flex.double([0.0] * n_bins * n_unique_indices)
    numerator.reshape(flex.grid(n_bins, n_unique_indices))
    denominator = flex.double([0.0] * n_unique_indices)
    for n in range(n_unique_indices):
        a1 = 0.0
        b1 = 0.0
        lsum = data_manager.h_index_counter_array[n]
        for i in range(0,lsum):
            indexer = i + data_manager.h_index_cumulative_array[n]
            l = data_manager.sorted_reflections['l_bin_index'][indexer]
            h = data_manager.sorted_reflections['h_index'][indexer]
            #if l != 0 and h != 0:
            numerator[l,h] += ((intensities[indexer]/(variances[indexer])) 
                                   -(2.0*g_values[l]*Ih_array[h]/(variances[indexer])))
            denominator[h] += ((g_values[l]**2)/variances[indexer])    
    for h in range(n_unique_indices):
        #if denominator[h] != 0.0:
        for l in range(n_bins):
            dIh_array[l,h] = numerator[l,h] / denominator[h]
        #else:
        #    print "zero denominator"
    #print len(dIh_array)
    #print dIh_array.count(0)
    return dIh_array