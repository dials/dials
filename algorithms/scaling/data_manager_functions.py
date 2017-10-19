'''
Define a Data_Manager object used for calculating scaling factors
'''
import copy
from dials.array_family import flex
from cctbx import miller, crystal
import minimiser_functions as mf
from dials.util.options import flatten_experiments, flatten_reflections
import numpy as np
import cPickle as pickle
from target_function import *
from basis_functions import *
from Wilson_outlier_test import calculate_wilson_outliers
import scale_factor as SF
from reflection_weighting import *
from data_quality_assessment import R_meas, R_pim
import matplotlib.pyplot as plt

class Data_Manager(object):
  '''Data Manager takes a params parsestring containing the parsed
     integrated.pickle and integrated_experiments.json files'''
  def __init__(self, reflections, experiments, scaling_options):
    'General attributes relevant for all parameterisations'
    self.experiments = experiments
    self.reflection_table = reflections[0]
    self.initial_keys = [key for key in self.reflection_table.keys()]
    self.reflection_table['centric_flag'] = flex.bool([False] * len(self.reflection_table))
    self.reflection_table['x_value'] = self.reflection_table['xyzobs.px.value'].parts()[0]
    self.reflection_table['y_value'] = self.reflection_table['xyzobs.px.value'].parts()[1]
    self.reflection_table['z_value'] = self.reflection_table['xyzobs.px.value'].parts()[2]
    self.reflection_table['inverse_scale_factor'] = flex.double([1.0] * len(self.reflection_table))
    self.reflection_table['Ih_values'] = flex.double([0.0] * len(self.reflection_table))
    #self.reflection_table['wilson_outlier_flag'] = flex.bool([False] * len(self.reflection_table))
    '''z_max = max(self.reflection_table['z_value'])
    sel = self.reflection_table['z_value'] < z_max/8.0
    self.reflection_table = self.reflection_table.select(sel)'''
    self.sorted_by_miller_index = False
    self.sorted_reflections = None
    self.reflections_for_scaling = None
    self.filtered_reflections = None
    self.Ih_array = None
    self.h_index_counter_array = None
    self.h_index_cumulative_array = None
    self.scaling_options = scaling_options
    #set choice of intensity and filter out bad values of d, variance etc.
    #fltering methods could be replaced by use of reflection table flags?
    self.select_optimal_intensities()
    self.filter_and_sort_reflections()
    self.update_weights_for_scaling(self.sorted_reflections)

  def extract_reflections_for_scaling(self, reflection_table):
    sel = self.weights_for_scaling.get_weights() > 0.0
    self.reflections_for_scaling = reflection_table.select(sel)
    self.assign_h_index(self.reflections_for_scaling)
    self.update_weights_for_scaling(self.reflections_for_scaling)

  def update_weights_for_scaling(self, reflection_table):
    self.weights_for_scaling = Weighting(reflection_table)
    self.weights_for_scaling.apply_Isigma_cutoff(reflection_table,
                                                 self.scaling_options['Isigma_min'])
    self.weights_for_scaling.apply_dmin_cutoff(reflection_table,
                                               self.scaling_options['d_min'])
    self.weights_for_scaling.remove_wilson_outliers(reflection_table)

  def filter_and_sort_reflections(self):
    self.filtered_reflections = copy.deepcopy(self.reflection_table)
    self.filter_negative_values(self.filtered_reflections, 'd')
    self.filter_negative_values(self.filtered_reflections, 'variance')
    self.filter_very_negative_intensities()#this may not be needed if we do some outlier rejection.
    self.map_indices_to_asu()
    self.assign_h_index(self.sorted_reflections)

  '''define a number of data filtering options'''
  def filter_negative_values(self, reflection_table, column):
    sel = reflection_table[column] > 0.0
    reflection_table = reflection_table.select(sel)

  def filter_very_negative_intensities(self):
    '''return boolean selection of a given variable range'''
    Ioversigma = self.filtered_reflections['intensity']/(self.filtered_reflections['variance']**0.5)
    sel = Ioversigma > -2.0
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

  def select_optimal_intensities(self):
    if (self.scaling_options['integration_method'] == 'sum' or
        self.scaling_options['integration_method'] == 'prf'):
      intstr = self.scaling_options['integration_method']
      self.reflection_table['intensity'] = (self.reflection_table['intensity.'+intstr+'.value']
                                            * self.reflection_table['lp']
                                            / self.reflection_table['dqe'])
      self.reflection_table['variance'] = (self.reflection_table['intensity.'+intstr+'.variance']
                                           * (self.reflection_table['lp']**2)
                                           / (self.reflection_table['dqe']**2))
    #option to add in future a combined prf/sum in a similar fashion to aimless
    elif self.scaling_options['integration_method'] == 'combine':
      self.scaling_options['integration_method'] = 'prf'
      self.select_optimal_intensities()

  def assign_h_index(self, reflection_table):
    '''assign an index to the sorted reflection table that
       labels each group of unique miller indices'''
    s = len(reflection_table['d'])
    if self.sorted_by_miller_index is False:
      raise ValueError('Data not yet sorted by miller index')
    else:
      reflection_table['h_index'] = flex.int([0] * s)
      self.h_index_counter_array = []
      h_index = 0
      h_index_counter = 1
      for i in range(1, s):
        if (reflection_table['asu_miller_index'][i] ==
            reflection_table['asu_miller_index'][i-1]):
          reflection_table['h_index'][i] = h_index
          h_index_counter += 1
        else:
          h_index += 1
          reflection_table['h_index'][i] = h_index
          self.h_index_counter_array.append(h_index_counter)
          h_index_counter = 1
      self.h_index_counter_array.append(h_index_counter)
      '''calculate the cumulative sum after each h_index group'''
      hsum = 0
      self.h_index_cumulative_array = [0]
      for n in self.h_index_counter_array:
        hsum += n
        self.h_index_cumulative_array.append(hsum)

  '''define a few methods for saving the data'''
  def save_sorted_reflections(self, filename):
    ''' Save the reflections to file. '''
    self.sorted_reflections.as_pickle(filename)

  def save_data_manager(self, filename):
    ''' Save the data manager to file. '''
    data_file = open(filename, 'w')
    pickle.dump(self, data_file)
    data_file.close()

  def calc_Ih(self, reflection_table):
    '''calculate the current best estimate for I for each reflection group'''
    intensities = reflection_table['intensity']
    variances = reflection_table['variance']
    scale_factors = reflection_table['inverse_scale_factor']
    scaleweights = self.weights_for_scaling.get_weights()
    #gsq = (((scale_factors)**2) / variances)
    gsq = (((scale_factors)**2) *scaleweights)
    sumgsq = flex.double(np.add.reduceat(gsq, self.h_index_cumulative_array[:-1]))
    #gI = ((scale_factors * intensities) / variances)
    gI = ((scale_factors * intensities) * scaleweights)
    sumgI = flex.double(np.add.reduceat(gI, self.h_index_cumulative_array[:-1]))
    sumweights = flex.double(np.add.reduceat(scaleweights, self.h_index_cumulative_array[:-1]))
    #self.Ih_array = sumgI / sumgsq
    self.Ih_array = flex.double([val/ sumgsq[i] if sumweights[i] > 0.0
                                 else 0.0 for i, val in enumerate(sumgI)])
    reflection_table['Ih_values'] = flex.double(
      np.repeat(self.Ih_array, self.h_index_counter_array))

  


class aimless_Data_Manager(Data_Manager):
  '''Data Manager subclass for implementing XDS parameterisation'''
  def __init__(self, reflections, experiments, scaling_options):
    Data_Manager.__init__(self, reflections, experiments, scaling_options)
    'Attributes specific to aimless parameterisation'
    '''set bin parameters'''
    self.binning_parameters = {'n_scale_bins' : None, 'n_B_bins' : None, 'lmax' : 4,
                               'n_d_bins' : None}
    self.bin_boundaries = None
    for key, value in scaling_options.iteritems():
      if key in self.binning_parameters:
        self.binning_parameters[key] = value
    #initialise g-value arrays
    self.g_absorption = None
    self.g_scale = None
    self.g_decay = None
    self.g_parameterisation = None
    self.active_parameterisation = None
    self.active_bin_index = None
    self.n_active_params = 0
    self.initialise_scale_factors()
    self.sorted_reflections['Wilson_outlier_flag'] = calculate_wilson_outliers(
      self.sorted_reflections, self.experiments)
    self.update_weights_for_scaling(self.sorted_reflections)
    self.reflections_for_scaling = copy.deepcopy(self.sorted_reflections)
    self.extract_reflections_for_scaling(self.sorted_reflections)
    #refactor the next two operations into extract_reflections?
    self.g_scale.set_normalised_values(self.reflections_for_scaling['normalised_rotation_angle'])
    self.g_decay.set_normalised_values(self.reflections_for_scaling['normalised_time_values'])


  def initialise_scale_factors(self):
    self.bin_reflections_scale()
    self.bin_reflections_decay()
    #self.initialise_absorption_scales()
    self.active_parameters = flex.double([])
    self.active_parameters.extend(self.g_scale.get_scale_factors())
    self.active_parameters.extend(self.g_decay.get_scale_factors())
    #self.active_parameters.extend(self.g_absorption)

  def calculate_scale_factors(self):
    scale_factors = self.g_scale.get_scale_factors()
    sorted_refl_scale_factor = SF.SmoothScaleFactor_1D(1.0, self.n_g_scale_params)
    sorted_refl_scale_factor.set_scale_factors(scale_factors)
    sorted_refl_scale_factor.set_normalised_values(self.sorted_reflections['normalised_rotation_angle'])
    self.sorted_reflections['angular_scale_factor'] = sorted_refl_scale_factor.calculate_smooth_scales()
    b_factors = self.g_decay.get_scale_factors()
    sorted_refl_decay_factor = SF.SmoothScaleFactor_1D(0.0, self.n_g_decay_params)
    sorted_refl_decay_factor.set_scale_factors(b_factors)
    sorted_refl_decay_factor.set_normalised_values(self.sorted_reflections['normalised_time_values'])
    b_factors = sorted_refl_decay_factor.calculate_smooth_scales()
    B_factor_values = flex.double(np.exp(b_factors/(2.0 * (self.sorted_reflections['d']**2))))
    self.sorted_reflections['decay_factor'] = B_factor_values
    self.sorted_reflections['inverse_scale_factor'] = (self.sorted_reflections['angular_scale_factor']
                                                       * self.sorted_reflections['decay_factor'])
     
  def get_target_function(self):
    '''call the aimless target function method'''
    return target_function(self).return_targets()

  def get_basis_function(self, parameters):
    '''call the aimless basis function method'''
    return aimless_basis_function(self, parameters).return_basis()

  def update_for_minimisation(self, parameters):
    '''update the scale factors and Ih for the next iteration of minimisation'''
    basis_fn = self.get_basis_function(parameters)
    self.reflections_for_scaling['inverse_scale_factor'] = basis_fn[0]
    self.active_derivatives = basis_fn[1]
    self.calc_Ih(self.reflections_for_scaling)

  def bin_reflections_decay(self):
    '''bin reflections for decay correction'''
    #still need to asign a bin index for wilson stats analysis.
    ndbins = self.binning_parameters['n_d_bins']
    resmax = (1.0 / (min(self.filtered_reflections['d'])**2))
    resmin = (1.0 / (max(self.filtered_reflections['d'])**2))
    resolution_bins = ((flex.double(range(0, ndbins + 1))
                        * ((resmax - resmin) / ndbins))
                       + flex.double([resmin] * (ndbins + 1)))
    d_bins = (1.0 /(resolution_bins[::-1]**0.5))
    #add a small tolerance to make sure no rounding errors cause extreme
    #data to not be selected
    d_bins[0] = d_bins[0] - 0.00001
    d_bins[-1] = d_bins[-1] + 0.00001
    firstbin_index = flex.int([-1] * len(self.sorted_reflections['d']))
    for i in range(ndbins):
      selection = select_variables_in_range(
        self.sorted_reflections['d'], d_bins[i], d_bins[i+1])
      firstbin_index.set_selected(selection, i)
    self.sorted_reflections['res_bin_index'] = firstbin_index

    rotation_interval = 15.0
    osc_range = self.experiments.scans()[0].get_oscillation_range()
    one_oscillation_width = self.experiments.scans()[0].get_oscillation()[1]
    #extend range by 0.001 to make sure all datapoints are within the min/max range
    (initial_rotation_angle, final_rotation_angle) = (osc_range[0] - 0.001, osc_range[1] + 0.001)
    self.sorted_reflections['normalised_time_values'] = ((self.sorted_reflections['z_value'] 
      * one_oscillation_width) - initial_rotation_angle)/rotation_interval
    #define the highest and lowest gridpoints: go out one further than the max/min int values
    highest_parameter_value = int((max(self.sorted_reflections['normalised_time_values'])//1)+3)#was +2
    lowest_parameter_value = int((min(self.sorted_reflections['normalised_time_values'])//1)-2)#was -1
    n_decay_parameters =  highest_parameter_value - lowest_parameter_value + 1
    self.g_decay = SF.SmoothScaleFactor_1D(0.0, n_decay_parameters)
    #self.g_decay.set_normalised_values(normalised_time_values)
    self.n_active_params += n_decay_parameters
    self.n_g_decay_params = n_decay_parameters

  def bin_reflections_scale(self):
    '''calculate the relative '''
    rotation_interval = 15.0
    osc_range = self.experiments.scans()[0].get_oscillation_range()
    one_oscillation_width = self.experiments.scans()[0].get_oscillation()[1]
    #extend range by 0.001 to make sure all datapoints are within the min/max range
    (initial_rotation_angle, final_rotation_angle) = (osc_range[0] - 0.001, osc_range[1] + 0.001)
    self.sorted_reflections['normalised_rotation_angle'] = ((self.sorted_reflections['z_value'] 
      * one_oscillation_width) - initial_rotation_angle)/rotation_interval
    #define the highest and lowest gridpoints: go out one further than the max/min int values
    highest_parameter_value = int((max(self.sorted_reflections['normalised_rotation_angle'])//1)+3)#was +2
    lowest_parameter_value = int((min(self.sorted_reflections['normalised_rotation_angle'])//1)-2)#was -1
    n_scale_parameters =  highest_parameter_value - lowest_parameter_value + 1
    #normalised_rotation_values = range(lowest_parameter_value, (highest_parameter_value+1))
    self.g_scale = SF.SmoothScaleFactor_1D(1.0, n_scale_parameters)
    #self.g_scale.set_normalised_values(normalised_rotation_values)
    self.n_active_params += n_scale_parameters
    self.n_g_scale_params = n_scale_parameters
    
  def initialise_absorption_scales(self):
    n_abs_params = 0
    for i in range(self.binning_parameters['lmax']):
      n_abs_params += (2*(i+1))+1
    self.g_absorption = SF.ScaleFactor(1.0, nabs_params)
    #self.g_absorption = flex.double([1.0] * n_abs_params)
    self.n_active_params += n_abs_params
    self.n_g_abs_params = n_abs_params

  def clean_reflection_table(self):
    #add keys for additional data that is to be exported
    self.initial_keys.append('inverse_scale_factor')
    for key in self.reflection_table.keys():
      if not key in self.initial_keys:
        del self.sorted_reflections[key]
    added_columns = ['Ih_values', 'h_index', 'asu_miller_index',
                     'res_bin_index', 'decay_factor', 'angular_scale_factor',
                     'normalised_rotation_angle', 'normalised_time_values']
    for key in added_columns:
      del self.sorted_reflections[key]

  


class XDS_Data_Manager(Data_Manager):
  '''Data Manager subclass for implementing XDS parameterisation'''
  def __init__(self, reflections, experiments, scaling_options):
    Data_Manager.__init__(self, reflections, experiments, scaling_options)
    'Attributes specific to XDS parameterisation'
    '''set bin parameters'''
    self.binning_parameters = {'n_d_bins' : None, 'n_z_bins' : None,
                               'n_absorption_positions' : 5,
                               'n_detector_bins' : None}
    self.bin_boundaries = None
    for key, value in scaling_options.iteritems():
      if key in self.binning_parameters:
        self.binning_parameters[key] = value
    #initialise g-value arrays
    self.g_absorption = None
    self.g_modulation = None
    self.g_decay = None
    self.g_parameterisation = None
    self.active_parameterisation = None
    self.active_bin_index = None
    self.n_active_params = None
    self.constant_g_values = None
    self.initialise_scale_factors()
    self.sorted_reflections['Wilson_outlier_flag'] = calculate_wilson_outliers(
      self.sorted_reflections, self.experiments)
    self.update_weights_for_scaling(self.sorted_reflections)
    self.reflections_for_scaling = copy.deepcopy(self.sorted_reflections)
    self.extract_reflections_for_scaling(self.sorted_reflections)
    self.g_modulation.set_normalised_values(self.reflections_for_scaling['normalised_x_values'],
      self.reflections_for_scaling['normalised_y_values'])
    self.g_decay.set_normalised_values(self.reflections_for_scaling['normalised_res_values'],
      self.reflections_for_scaling['normalised_time_values'])
    self.g_absorption.set_normalised_values(self.reflections_for_scaling['normalised_x_abs_values'],
      self.reflections_for_scaling['normalised_y_abs_values'], 
      self.reflections_for_scaling['normalised_time_values'])
    self.calculate_derivatives(self.reflections_for_scaling)  
    self.g_parameterisation = {
      'g_absorption' : {'index': 'a_bin_index', 'parameterisation' : self.g_absorption.get_scale_factors(),
                        'derivatives': self.g_absorption.calculate_smooth_derivatives(),
                        'values' : self.g_absorption.get_scales_of_reflections(),
                        'object' : self.g_absorption},
      'g_modulation' : {'index': 'xy_bin_index', 'parameterisation' : self.g_modulation.get_scale_factors(),
                        'derivatives': self.g_modulation.calculate_smooth_derivatives(),
                        'values' : self.g_modulation.get_scales_of_reflections(),
                        'object' : self.g_modulation},
      'g_decay' : {'index': 'l_bin_index', 'parameterisation' : self.g_decay.get_scale_factors(),
                   'derivatives': self.g_decay.calculate_smooth_derivatives(),
                   'values' : self.g_decay.get_scales_of_reflections(),
                   'object' : self.g_decay}}


  def calculate_derivatives(self, reflection_table):
    self.g_decay_derivatives = self.g_decay.calculate_smooth_derivatives()
    self.g_modulation_derivatives = self.g_modulation.calculate_smooth_derivatives()
    self.g_absorption_derivatives = self.g_absorption.calculate_smooth_derivatives()

  def initialise_scale_factors(self):
    self.bin_reflections_decay()
    self.bin_reflections_absorption()
    self.bin_reflections_modulation()

  def get_target_function(self):
    '''call the xds target function method'''
    if self.scaling_options['parameterization'] == 'log':
      return xds_target_function_log(self).return_targets()
    else:
      return target_function(self).return_targets()

  def get_basis_function(self, parameters):
    '''call the xds basis function method'''
    if self.scaling_options['parameterization'] == 'log':
      return xds_basis_function_log(self, parameters).return_basis()
    else:
      return xds_basis_function(self, parameters).return_basis()

  def update_for_minimisation(self, parameters):
    '''update the scale factors and Ih for the next iteration of minimisation'''
    basis_fn = self.get_basis_function(parameters)
    self.reflections_for_scaling['inverse_scale_factor'] = basis_fn[0]
    self.active_derivatives = basis_fn[1]
    self.calc_Ih(self.reflections_for_scaling)

  def calculate_scale_factors(self):
    modulation_factors = self.g_modulation.get_scale_factors()
    #print list(scale_factors)
    sorted_refl_mod_factor = SF.SmoothScaleFactor_2D(1.0, 
      self.g_modulation.n1_parameters, self.g_modulation.n2_parameters)
    sorted_refl_mod_factor.set_scale_factors(modulation_factors)
    sorted_refl_mod_factor.set_normalised_values(self.sorted_reflections['normalised_x_values'],
      self.sorted_reflections['normalised_y_values'])
    gscalevalues = sorted_refl_mod_factor.calculate_smooth_scales()

    decay_factors = self.g_decay.get_scale_factors()
    #print list(scale_factors)
    sorted_refl_decay_factor = SF.SmoothScaleFactor_2D(1.0, 
      self.g_decay.n1_parameters, self.g_decay.n2_parameters)
    sorted_refl_decay_factor.set_scale_factors(decay_factors)
    sorted_refl_decay_factor.set_normalised_values(self.sorted_reflections['normalised_res_values'],
      self.sorted_reflections['normalised_time_values'])
    gdecayvalues = sorted_refl_decay_factor.calculate_smooth_scales()

    abs_factors = self.g_absorption.get_scale_factors()
    #print list(scale_factors)
    sorted_refl_abs_factor = SF.SmoothScaleFactor_GridAbsorption(1.0, 
      self.g_absorption.nx_parameters, self.g_absorption.ny_parameters, 
      self.g_absorption.ntime_parameters)
    sorted_refl_abs_factor.set_scale_factors(abs_factors)
    sorted_refl_abs_factor.set_normalised_values(self.sorted_reflections['normalised_x_abs_values'],
      self.sorted_reflections['normalised_y_abs_values'], 
      self.sorted_reflections['normalised_time_values'])
    gabsvalues = sorted_refl_abs_factor.calculate_smooth_scales()

    if self.scaling_options['parameterization'] == 'log':
      self.sorted_reflections['inverse_scale_factor'] = flex.double(
        np.exp(gscalevalues + gdecayvalues + gabsvalues))
    else:
      self.sorted_reflections['inverse_scale_factor'] = (gscalevalues * gdecayvalues
                                                         * gabsvalues)

  def bin_reflections_decay(self):
    '''bin reflections for decay correction'''
    ndbins = self.binning_parameters['n_d_bins']
    nzbins = self.binning_parameters['n_z_bins']
    '''Bin the data into resolution and time 'z' bins'''
    zmax = max(self.filtered_reflections['z_value']) + 0.001
    zmin = min(self.filtered_reflections['z_value']) - 0.001

    resmax = (1.0 / (min(self.filtered_reflections['d'])**2)) + 0.001
    resmin = (1.0 / (max(self.filtered_reflections['d'])**2)) - 0.001

    one_resbin_width = (resmax - resmin) / ndbins
    one_time_width = (zmax - zmin) / nzbins

    self.sorted_reflections['normalised_res_values'] = (((1.0 / (self.sorted_reflections['d']**2)) 
      - resmin) / one_resbin_width)
    #define the highest and lowest gridpoints: go out one further than the max/min int values
    highest_parameter_value = int((max(self.sorted_reflections['normalised_res_values'])//1)+2)
    lowest_parameter_value = int((min(self.sorted_reflections['normalised_res_values'])//1)-1)
    n_res_parameters =  highest_parameter_value - lowest_parameter_value + 1
    #normalised_res_values = range(lowest_parameter_value, (highest_parameter_value+1))

    self.sorted_reflections['normalised_time_values'] = ((self.sorted_reflections['z_value'] 
      - zmin) / one_time_width)
    #define the highest and lowest gridpoints: go out one further than the max/min int values
    highest_parameter_value = int((max(self.sorted_reflections['normalised_time_values'])//1)+2)
    lowest_parameter_value = int((min(self.sorted_reflections['normalised_time_values'])//1)-1)
    n_time_parameters =  highest_parameter_value - lowest_parameter_value + 1
    #normalised_time_values = range(lowest_parameter_value, (highest_parameter_value+1))

    if self.scaling_options['parameterization'] == 'log':
      self.g_decay = SF.SmoothScaleFactor_2D(0.0, n_res_parameters, n_time_parameters)
    else:
      self.g_decay = SF.SmoothScaleFactor_2D(1.0, n_res_parameters, n_time_parameters)
    self.g_decay.set_normalised_values(self.sorted_reflections['normalised_res_values'], 
      self.sorted_reflections['normalised_time_values'])
  
    resolution_bins = ((flex.double(range(0, ndbins + 1))
                        * ((resmax - resmin) / ndbins))
                       + flex.double([resmin] * (ndbins + 1)))
    d_bins = (1.0 /(resolution_bins[::-1]**0.5))
    z_bins = ((flex.double(range(0, nzbins + 1)) * ((zmax - zmin) / nzbins))
              + flex.double([zmin] * (nzbins + 1)))
    #add a small tolerance to make sure no rounding errors cause extreme
    #data to not be selected
    d_bins[0] = d_bins[0] - 0.00001
    d_bins[-1] = d_bins[-1] + 0.00001
    z_bins[0] = z_bins[0] - 0.00001
    z_bins[-1] = z_bins[-1] + 0.00001
    self.bin_boundaries = {'d' : d_bins, 'z_value' : z_bins}

    firstbin_index = flex.int([-1] * len(self.sorted_reflections['d']))
    for i in range(ndbins):
      selection = select_variables_in_range(
        self.sorted_reflections['d'], d_bins[i], d_bins[i+1])
      firstbin_index.set_selected(selection, i)
    self.sorted_reflections['res_bin_index'] = firstbin_index


  def bin_reflections_modulation(self):
    '''bin reflections for modulation correction'''
    nxbins = nybins = self.binning_parameters['n_detector_bins']
    xvalues = self.sorted_reflections['x_value']
    (xmax, xmin) = (max(xvalues) + 0.001, min(xvalues) - 0.001)
    yvalues = self.sorted_reflections['y_value']
    (ymax, ymin) = (max(yvalues) + 0.001, min(yvalues) - 0.001)

    one_x_width = (xmax - xmin) / float(nxbins)
    one_y_width = (ymax - ymin) / float(nybins)

    self.sorted_reflections['normalised_x_values'] = ((xvalues - xmin) / one_x_width)
    self.sorted_reflections['normalised_y_values'] = ((yvalues - ymin) / one_y_width)
    #define the highest and lowest gridpoints: go out one further than the max/min int values
    highest_x_parameter_value = int((max(self.sorted_reflections['normalised_x_values'])//1)+2)
    lowest_x_parameter_value = int((min(self.sorted_reflections['normalised_x_values'])//1)-1)
    n_x_parameters =  highest_x_parameter_value - lowest_x_parameter_value + 1
    highest_y_parameter_value = int((max(self.sorted_reflections['normalised_y_values'])//1)+2)
    lowest_y_parameter_value = int((min(self.sorted_reflections['normalised_y_values'])//1)-1)
    n_y_parameters =  highest_y_parameter_value - lowest_y_parameter_value + 1

    if self.scaling_options['parameterization'] == 'log':
      #self.g_modulation = flex.double([0.0] * (nxbins ** 2))
      self.g_modulation = SF.SmoothScaleFactor_2D(0.0, n_x_parameters, n_y_parameters)
    else:
      #self.g_modulation = flex.double([1.0] * (nxbins ** 2))
      self.g_modulation = SF.SmoothScaleFactor_2D(1.0, n_x_parameters, n_y_parameters)
    self.g_modulation.set_normalised_values(self.sorted_reflections['normalised_x_values'], 
      self.sorted_reflections['normalised_y_values'])

  def bin_reflections_absorption(self):
    '''bin reflections for absorption correction'''
    nxbins = nybins = 4.0#self.binning_parameters['n_absorption_positions']
    nzbins = self.binning_parameters['n_z_bins']

    xvalues = self.sorted_reflections['x_value']
    (xmax, xmin) = (max(xvalues) + 0.001, min(xvalues) - 0.001)
    yvalues = self.sorted_reflections['y_value']
    (ymax, ymin) = (max(yvalues) + 0.001, min(yvalues) - 0.001)

    one_x_width = (xmax - xmin) / float(nxbins)
    one_y_width = (ymax - ymin) / float(nybins)

    self.sorted_reflections['normalised_x_abs_values'] = ((xvalues - xmin) / one_x_width)
    self.sorted_reflections['normalised_y_abs_values'] = ((yvalues - ymin) / one_y_width)
    #define the highest and lowest gridpoints: go out one further than the max/min int values
    highest_x_parameter_value = int((max(self.sorted_reflections['normalised_x_abs_values'])//1)+1)
    lowest_x_parameter_value = int((min(self.sorted_reflections['normalised_x_abs_values'])//1))
    n_x_parameters =  highest_x_parameter_value - lowest_x_parameter_value + 1
    highest_y_parameter_value = int((max(self.sorted_reflections['normalised_y_abs_values'])//1)+1)
    lowest_y_parameter_value = int((min(self.sorted_reflections['normalised_y_abs_values'])//1))
    n_y_parameters =  highest_y_parameter_value - lowest_y_parameter_value + 1

    n_time_bins = int((max(self.sorted_reflections['normalised_time_values'])//1)+1)

    if self.scaling_options['parameterization'] == 'log':
      self.g_absorption = SF.SmoothScaleFactor_GridAbsorption(0.0, n_x_parameters, n_y_parameters, n_time_bins)
    else:
      self.g_absorption = SF.SmoothScaleFactor_GridAbsorption(1.0, n_x_parameters, n_y_parameters, n_time_bins)
    self.g_absorption.set_normalised_values(self.sorted_reflections['normalised_x_abs_values'], 
      self.sorted_reflections['normalised_y_abs_values'], self.sorted_reflections['normalised_time_values'])


  def bin_reflections_absorption_radially(self):
    '''bin reflections for absorption correction'''
    from math import pi
    nxbins = nybins = self.binning_parameters['n_absorption_positions']
    nzbins = self.binning_parameters['n_z_bins']
    '''Bin the data into detector position and time 'z' bins'''
    z_bins = self.bin_boundaries['z_value']
    #define simple detector area map#
    image_size = self.experiments.detectors()[0].to_dict()['panels'][0]['image_size']
    xvalues = self.sorted_reflections['x_value']
    yvalues = self.sorted_reflections['y_value']
    x_center = image_size[0]/2.0
    y_center = image_size[1]/2.0
    radial_divider = max(x_center, y_center)
    xrelvalues = xvalues - x_center #!may need better definition of centerpoint
    yrelvalues = yvalues - y_center #!may need better definition of centerpoint
    radial_bins = [-0.0001, radial_divider / 3.0, 2.0 * radial_divider / 3.0,
                   ((((image_size[0]**2) + (image_size[1]**2))**0.5) / 2.0) + 5.0]
                   # '''+5.0 to the last bin adds extra tolerance to catch any
                   # spots with centres outside the detector area'''
    angular_bins = [-0.0001, pi / 4.0, 2.0 * pi / 4.0, 3.0 * pi / 4.0, pi,
                    5.0 * pi / 4.0, 6.0 * pi / 4.0, 7.0 * pi / 4.0, 2.0 * pi]
    radial_values = ((xrelvalues**2) + (yrelvalues**2))**0.5
    angular_values = flex.double(np.arccos(yrelvalues/radial_values))
    for i in range(0, len(angular_values)):
      if xrelvalues[i] < 0.0:
        angular_values[i] = (2.0 * pi) - angular_values[i]

    firstbin_index = flex.int([-1] * len(self.sorted_reflections['z_value']))
    secondbin_index = flex.int([-1] * len(self.sorted_reflections['z_value']))
    for i in range(len(angular_bins) - 1):
      selection1 = select_variables_in_range(angular_values, angular_bins[i],
                                             angular_bins[i+1])
      for j in range(len(radial_bins) - 1):
        selection2 = select_variables_in_range(radial_values, radial_bins[j],
                                               radial_bins[j+1])
        firstbin_index.set_selected(selection1 & selection2,
                                    ((i * (len(radial_bins) - 1)) + j))
    for i in range(nzbins):
      selection = select_variables_in_range(self.sorted_reflections['z_value'],
                                            z_bins[i], z_bins[i+1])
      secondbin_index.set_selected(selection, i)

    if firstbin_index.count(-1) > 0 or secondbin_index.count(-1) > 0:
      raise ValueError('Unable to fully bin data for absorption in scaling initialisation')
    self.sorted_reflections['a_bin_index'] = (firstbin_index
      + (secondbin_index * (len(angular_bins) - 1) * (len(radial_bins) - 1)))
    self.n_absorption_bins = (len(angular_bins) - 1) * (len(radial_bins) - 1)
    if self.scaling_options['parameterization'] == 'log':
      #self.g_absorption = flex.double([0.0] * self.n_absorption_bins * nzbins)
      self.g_absorption = SF.ScaleFactor(0.0, (self.n_absorption_bins * nzbins))
    else:
      #self.g_absorption = flex.double([1.0] * self.n_absorption_bins * nzbins)
      self.g_absorption = SF.ScaleFactor(1.0, (self.n_absorption_bins * nzbins))


  def scale_gvalues(self):
    '''Rescale the decay g-values by a relative B-factor and a global scale
    factor. '''
    Optimal_rescale_values = mf.B_optimiser(self, flex.double([0.0, 1.0]))
    print "Brel, 1/global scale = "+str(list(Optimal_rescale_values.x))

    scaling_factors = []
    for _ in range(self.binning_parameters['n_z_bins']):
      scaling_factors += flex.exp(Optimal_rescale_values.x[0]
                    *Optimal_rescale_values.res_values)
    scaling_factors = flex.double(scaling_factors)

    self.g_decay = self.g_decay * scaling_factors
    self.g_decay = self.g_decay * (1.0 / Optimal_rescale_values.x[1])
    print "scaled by B_rel and global scale parameter"

  def create_fake_dataset(self):
    self.sorted_reflections['intensity.sum.value'] = (
      self.sorted_reflections['Ih_values'] *
      self.sorted_reflections['inverse_scale_factor'] *
      self.sorted_reflections['dqe'] / self.sorted_reflections['lp'])

  def clean_reflection_table(self):
    #add keys for additional data that is to be exported
    self.initial_keys.append('inverse_scale_factor')
    for key in self.reflection_table.keys():
      if not key in self.initial_keys:
        del self.sorted_reflections[key]
    added_columns = ['l_bin_index', 'a_bin_index', 'xy_bin_index', 'h_index',
                     'asu_miller_index', 'res_bin_index', 'normalised_y_values',
                     'normalised_x_values', 'normalised_y_abs_values',
                     'normalised_x_abs_values', 'normalised_time_values', 
                     'normalised_res_values']
    for key in added_columns:
      del self.sorted_reflections[key]





def select_variables_in_range(variable_array, lower_limit, upper_limit):
  '''return boolean selection of a given variable range'''
  sel = flex.bool()
  for variable in variable_array:
    if lower_limit < variable <= upper_limit:
      sel.append(True)
    else:
      sel.append(False)
  return sel