'''
Define a Data_Manager object used for calculating scaling factors
'''
import copy
from dials_array_family_flex_ext import *
from cctbx.array_family.flex import *
from cctbx.array_family import flex
from cctbx import miller, crystal
import minimiser_functions as mf
from dials.util.options import (flatten_experiments, flatten_reflections)
import numpy as np

class Data_Manager(object):
    '''Data Manager takes a params parsestring
       containing the parsed integrated.pickle
       and integrated_experiments.json files'''
    def __init__(self, params, int_method):
        'General attributes relevant for all parameterisations'
        self.experiments = flatten_experiments(params.input.experiments)
        self.reflection_table = flatten_reflections(params.input.reflections)[0]
        n_entries = range(len(self.reflection_table['xyzobs.px.value']))
        self.reflection_table['x_value'] = flex.double(
            [self.reflection_table['xyzobs.px.value'][i][0] for i in n_entries])
        self.reflection_table['y_value'] = flex.double(
            [self.reflection_table['xyzobs.px.value'][i][1] for i in n_entries])
        self.reflection_table['z_value'] = flex.double(
            [self.reflection_table['xyzobs.px.value'][i][2] for i in n_entries])
        self.reflection_table['resolution'] = flex.double(
            [1.0/(x**2) if x != 0 else 0 for x in self.reflection_table['d']])
        self.filtered_reflections = copy.deepcopy(self.reflection_table)
        self.int_method = int_method
        self.sorted_by_miller_index = False
        self.sorted_reflections = None
        self.Ih_array = None
        self.Ih_values = None
        self.h_index_counter_array = None
        self.h_index_cumulative_array = None

    def filter_data(self, reflection_table_key, lower, upper):
        '''Filter reflection data for a given measurement variable and limits'''
        bad_data = select_variables_in_range(
            self.filtered_reflections[reflection_table_key], lower, upper)
        inv_sel = ~bad_data
        self.filtered_reflections = self.filtered_reflections.select(inv_sel)

    def filter_I_sigma(self, ratio):
        sel = flex.bool()
        for index, intensity in enumerate(self.filtered_reflections[self.int_method[0]]):
            if intensity/(self.filtered_reflections[self.int_method[1]][index]**0.5) > ratio:
                sel.append(True)
            else:
                sel.append(False)
        self.filtered_reflections = self.filtered_reflections.select(sel)


    def map_indices_to_asu(self):
        '''Create a miller_set object, map to the asu and create a sorted
           reflection table, sorted by asu miller index'''
        u_c = self.experiments.crystals()[0].get_unit_cell().parameters()
        s_g = self.experiments.crystals()[0].get_space_group()
        crystal_symmetry = crystal.symmetry(unit_cell=u_c, space_group=s_g)
        miller_set = miller.set(crystal_symmetry=crystal_symmetry,
                                indices=self.filtered_reflections['miller_index'])
        self.filtered_reflections["asu_miller_index"] = miller_set.map_to_asu().indices()
        permuted = (miller_set.map_to_asu()).sort_permutation(by_value='packed_indices')
        self.sorted_reflections = self.filtered_reflections.select(permuted)
        self.sorted_by_miller_index = True

    def scale_by_LP_and_dqe(self):
        '''Apply Lorenz polarisation and dqe correction to intensities
        and variances'''
        for q in self.int_method:
            if q.split('.')[2] == 'value':
                self.sorted_reflections[q+'.LPscaled'] = (self.sorted_reflections[q]
                                                          * self.sorted_reflections['lp']
                                                          / self.sorted_reflections['dqe'])
            elif q.split('.')[2] == 'variance':
                self.sorted_reflections[q+'.LPscaled'] = (self.sorted_reflections[q]
                                                          * (self.sorted_reflections['lp']**2)
                                                          / (self.sorted_reflections['dqe']**2))
        new_int_method = ['.', '.']
        for i, val in enumerate(self.int_method):
            new_int_method[i] = val+'.LPscaled'
        self.int_method = new_int_method

    def assign_h_index(self):
        '''assign an index to the sorted reflection table that
           labels each group of unique miller indices'''
        s = len(self.sorted_reflections['d'])
        if self.sorted_by_miller_index is False:
            raise ValueError('Data not yet sorted by miller index')
        else:
            self.sorted_reflections['h_index'] = flex.int([0] * s)
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

    def save_sorted_reflections(self, filename):
        ''' Save the reflections to file. '''
        self.sorted_reflections.as_pickle(filename)


class Kabsch_Data_Manager(Data_Manager):
    '''Data Manager subclass for implementing Kabsch parameterisation'''
    def __init__(self, params, int_method, gridding_parameters):
        Data_Manager.__init__(self, params, int_method)
        'Attributes specific to kabsch parameterisation'
        '''set bin parameters'''
        self.binning_parameters = {'ndbins':None, 'nzbins':None,
                                   'n_absorption_positions':5,
                                   'n_detector_bins':None}
        self.bin_boundaries = None
        for key, value in gridding_parameters.iteritems():
            self.binning_parameters[key] = value
        '''initialise g parameters'''
        self.scale_factors = None
        self.g_values = None
        self.g2_values = None
        self.g3_values = None

    def bin_reflections_decay(self):
        '''bin reflections for decay correction'''
        ndbins = self.binning_parameters['ndbins']
        nzbins = self.binning_parameters['nzbins']
        '''Bin the data into resolution and time 'z' bins'''
        zmax = max(self.filtered_reflections['z_value'])
        zmin = min(self.filtered_reflections['z_value'])
        resmax = (1.0 / (min(self.filtered_reflections['d'])**2))
        resmin = (1.0 / (max(self.filtered_reflections['d'])**2))
        resolution_bins = ((flex.double(range(0, ndbins + 1))
                            * ((resmax - resmin)/ndbins))
                           + flex.double([resmin] * (ndbins + 1)))
        d_bins = (1.0/(resolution_bins[::-1]**0.5))
        z_bins = ((flex.double(range(0, nzbins + 1)) * ((zmax - zmin)/nzbins))
                  + flex.double([zmin] * (nzbins + 1)))
        z_bins[0] = z_bins[0]-0.001
        firstbin_index = flex.int([-1] * len(self.sorted_reflections['d']))
        secondbin_index = flex.int([-1] * len(self.sorted_reflections['z_value']))
        for i in range(ndbins):
            selection = select_variables_in_range(
                self.sorted_reflections['d'], d_bins[i], d_bins[i+1])
            firstbin_index.set_selected(selection, i)
        for i in range(nzbins):
            selection = select_variables_in_range(
                self.sorted_reflections['z_value'], z_bins[i], z_bins[i+1])
            secondbin_index.set_selected(selection, i)
        self.sorted_reflections['l_bin_index'] = (firstbin_index
                                                  + (secondbin_index * ndbins))
        self.bin_boundaries = {'d' : d_bins, 'z_value' : z_bins}
        self.g_values = flex.double([1.0] * (ndbins * nzbins))

    def bin_reflections_modulation(self):
        '''bin reflections for modulation correction'''
        ngridpoints = self.binning_parameters['n_detector_bins']
        nxbins = nybins = ngridpoints
        xvalues = self.sorted_reflections['x_value']
        (xmax, xmin) = (max(xvalues), min(xvalues))
        yvalues = self.sorted_reflections['y_value']
        (ymax, ymin) = (max(yvalues), min(yvalues))
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
        self.sorted_reflections['xy_bin_index'] = (firstbin_index
                                                   + (secondbin_index * ngridpoints))
        self.g3_values = flex.double([1.0] * (ngridpoints ** 2))

    def bin_reflections_absorption(self):
        '''bin reflections for absorption correction'''
        nxbins = nybins = self.binning_parameters['n_absorption_positions']
        nzbins = self.binning_parameters['nzbins']
        '''Bin the data into detector position and time 'z' bins'''
        z_bins = self.bin_boundaries['z_value']
        #define simple detector area map#
        xvalues = self.sorted_reflections['x_value']
        (xmax, xmin) = (max(xvalues), min(xvalues))
        yvalues = self.sorted_reflections['y_value']
        (ymax, ymin) = (max(yvalues), min(yvalues))
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
        for i in range(nzbins):
            selection = select_variables_in_range(
                self.sorted_reflections['z_value'], z_bins[i], z_bins[i+1])
            secondbin_index.set_selected(selection, i)
        self.sorted_reflections['a_bin_index'] = (firstbin_index
                                                  + (secondbin_index * nxbins * nybins))
        self.g2_values = flex.double([1.0] * (nxbins * nybins * nzbins))

    def initialise_scale_factors(self):
        self.bin_reflections_decay()
        self.bin_reflections_absorption()
        self.bin_reflections_modulation()
        self.scale_factors = flex.double([1.0]*len(self.sorted_reflections['d']))

    def calc_Ih(self):
        '''calculate the current best estimate for I for each reflection group'''
        intensities = self.sorted_reflections[self.int_method[0]]
        variances = self.sorted_reflections[self.int_method[1]]
        gsq = (((self.scale_factors)**2)/variances)
        sumgsq = flex.double(np.add.reduceat(gsq, self.h_index_cumulative_array[:-1]))
        gI = ((self.scale_factors*intensities)/variances)
        sumgI = flex.double(np.add.reduceat(gI, self.h_index_cumulative_array[:-1]))
        self.Ih_array = sumgI/sumgsq
        self.Ih_values = flex.double(np.repeat(self.Ih_array, self.h_index_counter_array))

    def scale_gvalues(self):
        '''Rescale the decay g-values by a relative B-factor and a global scale
        factor. '''
        Optimal_rescale_values = mf.B_optimiser(self, flex.double([0.0, 1.0]))
        print "Brel, 1/global scale = "+str(list(Optimal_rescale_values.x))

        scaling_factors = []
        for _ in range(self.binning_parameters['nzbins']):
            scaling_factors += flex.exp(Optimal_rescale_values.x[0]
                                        *Optimal_rescale_values.res_values)
        scaling_factors = flex.double(scaling_factors)

        self.g_values = self.g_values * scaling_factors
        self.g_values = self.g_values * (1.0/Optimal_rescale_values.x[1])
        print "scaled by B_rel and global scale parameter"

    '''def clean_sorted_reflections(self):
        keylist_reflection_table = ['inverse_scale_factor']
        keylist_sorted_reflections = []
        self.sorted_reflections['inverse_scale_factor'] = self.scale_factors
        print list(self.sorted_reflections)[0]
        for key in self.reflection_table.keys():
            keylist_reflection_table.append(key)
        self.sorted_reflections = self.sorted_reflections[keylist_reflection_table]
        print list(self.sorted_reflections)[0]
        exit()'''


def select_variables_in_range(variable_array, lower_limit, upper_limit):
    '''return boolean selection of a given variable range'''
    sel = flex.bool()
    for variable in variable_array:
        if lower_limit < variable <= upper_limit:
            sel.append(True)
        else:
            sel.append(False)
    return sel
