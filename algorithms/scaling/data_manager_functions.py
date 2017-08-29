import copy
from dials_array_family_flex_ext import *
from cctbx.array_family.flex import *
from cctbx.array_family import flex
from cctbx import miller
from cctbx import crystal
import minimiser_functions as mf
from dials.util.options import (flatten_experiments, flatten_reflections)



class Data_Manager(object):
    '''Data Manager takes a params parsestring
       containing the parsed integrated.pickle
       and integrated_experiments.json files'''
    def __init__(self, params):
        self.experiments = flatten_experiments(params.input.experiments)
        self.reflection_table = flatten_reflections(params.input.reflections)[0]
        n_entries = len(self.reflection_table['xyzobs.px.value'])
        zvalues = flex.double([self.reflection_table['xyzobs.px.value'][i][2]
                               for i in range(n_entries)])
        self.reflection_table['z_value'] = zvalues
        self.reflection_table['resolution'] = flex.double(
            [1.0/(x**2) if x != 0 else 0 for x in self.reflection_table['d']])
        self.filtered_reflections = copy.deepcopy(self.reflection_table)
        self.sorted_by_miller_index = False
        self.sorted_reflections = None
        self.l_bin_index = None
        self.bin_boundaries = None
        self.g_values = None
        self.n_bins = None
        self.h_index_counter_array = None
        self.h_index_cumulative_array = None
        self.n_unique_indices = None
        self.Ih_array = None
        self.nzbins = None
        #repackage some of these attributes for conciseness

    def filter_data(self, reflection_table_key, lower, upper):
        '''Filter reflection data for a given measurement variable and limits'''
        bad_data = select_variables_in_range(self.filtered_reflections[reflection_table_key],
                                             lower, upper)
        inv_sel = ~bad_data
        self.filtered_reflections = self.filtered_reflections.select(inv_sel)

    def map_indices_to_asu(self):
        '''Create a miller_set object, map to the asu and create a sorted
           reflection table, sorted by asu miller index'''
        u_c = self.experiments.crystals()[0].get_unit_cell().parameters()
        s_g = self.experiments.crystals()[0].get_space_group()
        crystal_symmetry = crystal.symmetry(unit_cell=u_c, space_group=s_g)
        miller_set = miller.set(crystal_symmetry=crystal_symmetry,
                                indices=self.filtered_reflections['miller_index'])
        miller_array = miller.array(miller_set)
        self.filtered_reflections["asu_miller_index"] = miller_set.map_to_asu().indices()
        permuted = (miller_set.map_to_asu()).sort_permutation(by_value='packed_indices')
        self.sorted_reflections = self.filtered_reflections.select(permuted)
        self.sorted_by_miller_index = True

    def bin_reflections(self, bin_parameters):
        '''Takes a bin_parameters dict consisting of
        [[[data to bin],[list of bin boundaries]],....]'''
        firstbin_index = flex.int([0]*len(self.sorted_reflections[bin_parameters[0][0]]))
        secondbin_index = flex.int([0]*len(self.sorted_reflections[bin_parameters[1][0]]))
        for i in range(len(bin_parameters[0][1])-1):
            selection = select_variables_in_range(
                self.sorted_reflections[bin_parameters[0][0]],
                bin_parameters[0][1][i], bin_parameters[0][1][i+1])
            firstbin_index.set_selected(selection, i)
        for i in range(len(bin_parameters[1][1])-1):
            selection = select_variables_in_range(
                self.sorted_reflections[bin_parameters[1][0]],
                bin_parameters[1][1][i], bin_parameters[1][1][i+1])
            secondbin_index.set_selected(selection, i)
        self.l_bin_index = (firstbin_index
                            + (((secondbin_index))
                               * (len(bin_parameters[0][1])-1)))
        self.sorted_reflections['l_bin_index'] = self.l_bin_index
        self.bin_boundaries = {bin_parameters[0][0] : bin_parameters[0][1],
                               bin_parameters[1][0] : bin_parameters[1][1]}
        self.nzbins = len(bin_parameters[1][1]) - 1
        self.ndbins = len(bin_parameters[0][1]) - 1

    def set_g_values(self, gvalues):
        self.g_values = gvalues
        self.n_bins = len(gvalues)

    def assign_h_index(self):
        '''assign an index to the sorted reflection table that
           labels each group of unique miller indices'''
        s = len(self.sorted_reflections['d'])
        if self.sorted_by_miller_index is False:
            raise ValueError('Data not yet sorted by miller index')
        else:
            self.sorted_reflections['h_index'] = flex.int([0]*s)
            self.h_index_counter_array = []
            h_index = 0
            h_index_counter = 1
            for i in range(1, s):
                if  (self.sorted_reflections['asu_miller_index'][i] ==
                     self.sorted_reflections['asu_miller_index'][i-1]):
                    self.sorted_reflections['h_index'][i] = h_index
                    h_index_counter += 1
                else:
                    h_index += 1
                    self.sorted_reflections['h_index'][i] = h_index
                    self.h_index_counter_array.append(h_index_counter)
                    h_index_counter = 1
            self.h_index_counter_array.append(h_index_counter)
            '''calculate the cumulative sum after each h_index group'''
            hsum = 0
            self.h_index_cumulative_array = [0]
            for n in self.h_index_counter_array:
                hsum += n
                self.h_index_cumulative_array.append(hsum)
            self.n_unique_indices = len(self.h_index_counter_array)

    def calc_Ih(self):
        intensities = self.sorted_reflections['intensity.sum.value']
        g_values = self.g_values
        variances = self.sorted_reflections['intensity.sum.variance']
        self.Ih_array = []
        for h in range(self.n_unique_indices):
            a1 = 0.0
            b1 = 0.0
            lsum = self.h_index_counter_array[h]
            for i in range(lsum):
                indexer = i + self.h_index_cumulative_array[h]
                l = self.sorted_reflections['l_bin_index'][indexer]
                a1 += (g_values[l]*intensities[indexer]/variances[indexer])
                b1 += ((g_values[l]**2)/variances[indexer])
            self.Ih_array.append(a1/b1)

    def scale_gvalues(self):
        Optimal_rescale_values = mf.B_optimiser(self, flex.double([0.7, 10.0]))
        print list(Optimal_rescale_values.x)

        scaling_factors = []
        for i in range(0, self.nzbins):
            scaling_factors += flex.exp(Optimal_rescale_values.x[0]*Optimal_rescale_values.res_values)
        scaling_factors = flex.double(scaling_factors)

        self.g_values = self.g_values * scaling_factors
        self.g_values = self.g_values * (1.0/Optimal_rescale_values.x[1])
        self.calc_Ih()
        print "scaled by Brel and global scale parameter"


def select_variables_in_range(variable_array, lower_limit, upper_limit):
    sel = flex.bool()
    for j in range(len(variable_array)):
        if lower_limit < variable_array[j] <= upper_limit:
            sel.append(True)
        else:
            sel.append(False)
    return sel

def calculate_evenly_populated_resolution_bins(Loaded_reflections, ndbins):
    permuted = Loaded_reflections.millerset.map_to_asu().sort_permutation(
        by_value='resolution', reverse=True)
    print list(permuted)
    number_of_entries = len(Loaded_reflections.filtered_reflections['d'])
    sorted_by_d = (Loaded_reflections.filtered_reflections).select(permuted)
    print list(sorted_by_d['d'])[0:10]
    resolution_bins = flex.double([0]*(ndbins+1))
    for i in range(ndbins):
        n = (i*number_of_entries//ndbins)
        resolution_bins[i] = sorted_by_d['d'][n]
    resolution_bins[ndbins] = sorted_by_d['d'][-1]
    return resolution_bins
