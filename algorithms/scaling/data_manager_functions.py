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
    def __init__(self, params, int_str, var_str):
        self.experiments = flatten_experiments(params.input.experiments)
        self.reflection_table = flatten_reflections(params.input.reflections)[0]
        n_entries = len(self.reflection_table['xyzobs.px.value'])
        xvalues = flex.double([self.reflection_table['xyzobs.px.value'][i][0]
                               for i in range(n_entries)])
        yvalues = flex.double([self.reflection_table['xyzobs.px.value'][i][1]
                               for i in range(n_entries)])
        zvalues = flex.double([self.reflection_table['xyzobs.px.value'][i][2]
                               for i in range(n_entries)])
        self.reflection_table['x_value'] = xvalues
        self.reflection_table['y_value'] = yvalues
        self.reflection_table['z_value'] = zvalues
        self.reflection_table['resolution'] = flex.double(
            [1.0/(x**2) if x != 0 else 0 for x in self.reflection_table['d']])
        self.filtered_reflections = copy.deepcopy(self.reflection_table)
        self.int_str = int_str
        self.var_str = var_str
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
        bad_data = select_variables_in_range(
                   self.filtered_reflections[reflection_table_key],lower, upper)
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

    def scale_by_LP_and_dqe(self):
        self.sorted_reflections[self.int_str] = (self.sorted_reflections[self.int_str]
            * self.sorted_reflections['lp'] * self.sorted_reflections['dqe'])
        self.sorted_reflections[self.var_str] = (self.sorted_reflections[self.var_str]
            * self.sorted_reflections['lp'] * self.sorted_reflections['dqe'])

    def bin_reflections_dz(self, ndbins, nzbins):
        '''Bin the data into resolution and time 'z' bins'''
        zmax = max(self.filtered_reflections['z_value'])
        zmin = min(self.filtered_reflections['z_value'])
        resmax = (1.0 / (min(self.filtered_reflections['d'])**2))
        resmin = (1.0 / (max(self.filtered_reflections['d'])**2))
        resolution_bins = ((flex.double(range(0, ndbins + 1))
                            *((resmax - resmin)/ndbins))
                           +flex.double([resmin] * (ndbins + 1)))
        d_bins = (1.0/(resolution_bins[::-1]**0.5))
        z_bins = flex.double(range(0, nzbins + 1)) * zmax/nzbins
        firstbin_index = flex.int([-1]*len(self.sorted_reflections['d']))
        secondbin_index = flex.int([-1]*len(self.sorted_reflections['z_value']))
        for i in range(ndbins):
            selection = select_variables_in_range(
                self.sorted_reflections['d'], d_bins[i], d_bins[i+1])
            firstbin_index.set_selected(selection, i)
        for i in range(nzbins):
            selection = select_variables_in_range(
                self.sorted_reflections['z_value'], z_bins[i], z_bins[i+1])
            secondbin_index.set_selected(selection, i)
        self.l_bin_index = firstbin_index + (secondbin_index * ndbins)
        self.sorted_reflections['l_bin_index'] = self.l_bin_index
        self.bin_boundaries = {'d' : d_bins, 'z_value' : z_bins}
        self.nzbins = nzbins
        self.ndbins = ndbins
        self.g_values = flex.double([1.0] * (ndbins * nzbins))

    def bin_reflections_modulation(self,ngridpoints):
        nxbins = nybins = ngridpoints
        xvalues = self.sorted_reflections['x_value']
        xmax = max(xvalues); xmin = min(xvalues)
        yvalues = self.sorted_reflections['y_value']
        ymax = max(yvalues); ymin = min(yvalues)
        x_bins = (((flex.double(range(0, nxbins + 1)) * (xmax - xmin)/(nxbins)))
                  + flex.double([xmin] * (nxbins + 1)))
        x_bins[0] = x_bins[0]-0.001
        y_bins = ((flex.double(range(0, nybins + 1)) * (ymax - ymin)/(nybins)) 
                  + flex.double([ymin] * (nybins + 1)))
        y_bins[0] = y_bins[0]-0.001
        firstbin_index = flex.int([-1]*len(self.sorted_reflections['z_value']))
        secondbin_index = flex.int([-1]*len(self.sorted_reflections['z_value']))
        for i in range(nxbins):
            selection = select_variables_in_range(
                self.sorted_reflections['x_value'], x_bins[i], x_bins[i+1]) 
            firstbin_index.set_selected(selection, i)
        for i in range(nybins):
            selection = select_variables_in_range(
                self.sorted_reflections['y_value'], y_bins[i], y_bins[i+1]) 
            secondbin_index.set_selected(selection, i)
        xy_bin_index = firstbin_index + (secondbin_index * ngridpoints)
        self.sorted_reflections['xy_bin_index'] = xy_bin_index
        self.g3_values = flex.double([1.0] * (ngridpoints ** 2))
        self.ngridpoints = ngridpoints

    def bin_reflections_absorption(self, npos):
        nxbins = nybins = npos
        '''Bin the data into detector position and time 'z' bins'''
        z_bins = self.bin_boundaries['z_value']
        #define simple detector area map#
        xvalues = self.sorted_reflections['x_value']
        xmax = max(xvalues); xmin = min(xvalues)
        yvalues = self.sorted_reflections['y_value']
        ymax = max(yvalues); ymin = min(yvalues)

        x_bins = ((flex.double(range(0, nxbins + 1)) * (xmax - xmin)/(nxbins))
                  + flex.double([xmin] * (nxbins + 1)))
        x_bins[0] = x_bins[0]-0.001
        y_bins = ((flex.double(range(0, nybins + 1)) * (ymax - ymin)/(nybins))
                  + flex.double([ymin] * (nybins + 1)))
        y_bins[0] = y_bins[0]-0.001
        firstbin_index = flex.int([-1]*len(self.sorted_reflections['z_value']))
        secondbin_index = flex.int([-1]*len(self.sorted_reflections['z_value']))
        for i in range(nxbins):
            selection1 = select_variables_in_range(
                self.sorted_reflections['x_value'], x_bins[i], x_bins[i+1])
            for j in range(nybins):
                selection2 = select_variables_in_range(
                    self.sorted_reflections['y_value'], y_bins[j], y_bins[j+1])
                firstbin_index.set_selected(selection1 & selection2,
                                            ((i*nybins) + j))
        for i in range(self.nzbins):
            selection = select_variables_in_range(
                self.sorted_reflections['z_value'], z_bins[i], z_bins[i+1])
            secondbin_index.set_selected(selection, i)
        a_bin_index = firstbin_index + (secondbin_index * nxbins * nybins)
        self.sorted_reflections['a_bin_index'] = a_bin_index
        self.g2_values = flex.double([1.0] * (nxbins * nybins * self.nzbins))
        self.npos = npos

    def create_index_tuple(self):
        h_array = self.sorted_reflections['h_index']
        l_array = self.sorted_reflections['l_bin_index']
        a_array = self.sorted_reflections['a_bin_index']
        xy_array = self.sorted_reflections['xy_bin_index']
        self.bin_index_array = zip(h_array, l_array, a_array, xy_array)

    def set_g_values(self, gvalues):
        self.g_values = gvalues

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
        intensities = self.sorted_reflections[self.int_str]
        variances = self.sorted_reflections[self.var_str]
        self.Ih_array = []
        for h in range(self.n_unique_indices):
            a1 = 0.0
            b1 = 0.0
            lsum = self.h_index_counter_array[h]
            for i in range(lsum):
                indexer = i + self.h_index_cumulative_array[h]
                l = self.sorted_reflections['l_bin_index'][indexer]
                a = self.sorted_reflections['a_bin_index'][indexer]
                xy = self.sorted_reflections['xy_bin_index'][indexer]
                a1 += (self.g_values[l]*self.g2_values[a]*self.g3_values[xy]
                       *intensities[indexer]/variances[indexer])
                b1 += (((self.g_values[l]*self.g2_values[a]*self.g3_values[xy])**2)
                       /variances[indexer])
            self.Ih_array.append(a1/b1)

    def scale_gvalues(self):
        Optimal_rescale_values = mf.B_optimiser(self, flex.double([0.0, 1.0]))
        print "Brel, 1/global scale = "+str(list(Optimal_rescale_values.x))

        scaling_factors = []
        for i in range(0, self.nzbins):
            scaling_factors += flex.exp(Optimal_rescale_values.x[0]
                                        *Optimal_rescale_values.res_values)
        scaling_factors = flex.double(scaling_factors)

        self.g_values = self.g_values * scaling_factors
        self.g_values = self.g_values * (1.0/Optimal_rescale_values.x[1])
        self.calc_Ih()
        print "scaled by B_rel and global scale parameter"


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
