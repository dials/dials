'''
Define a Data_Manager object used for calculating scaling factors
'''
from __future__ import print_function
import copy
from dials.array_family import flex
from cctbx import miller, crystal
import minimiser_functions as mf
from dials.util.options import flatten_experiments, flatten_reflections
import numpy as np
import cPickle as pickle
from target_function import target_function, target_function_fixedIh, xds_target_function_log
import basis_functions as bf
from scaling_utilities import calc_s2d, sph_harm_table
from Wilson_outlier_test import calculate_wilson_outliers, calc_normE2
import scale_factor as SF
from reflection_weighting import Weighting
from data_quality_assessment import R_meas, R_pim
from target_Ih import single_Ih_table, joined_Ih_table, base_Ih_table
import minimiser_functions as mf
from collections import OrderedDict
from scitbx import sparse


class Data_Manager(object):
  '''Data Manager takes a params parsestring containing the parsed
     integrated.pickle and integrated_experiments.json files'''
  def __init__(self, reflections, experiments, scaling_options):
    'General attributes relevant for all parameterisations'
    print('\nInitialising a data manager instance. \n')
    self.experiments = experiments
    'initial filter to select integrated reflections'
    reflections = reflections.select(reflections['intensity.prf.variance'] > 0)
    reflections = reflections.select(reflections['intensity.sum.variance'] > 0)
    self.reflection_table = reflections.select(reflections.get_flags(
      reflections.flags.integrated))
    self.scaling_options = scaling_options
    self.initial_keys = [key for key in self.reflection_table.keys()]
    if not 'inverse_scale_factor' in self.initial_keys:
      self.reflection_table['inverse_scale_factor'] = (
        flex.double([1.0] * len(self.reflection_table)))
    self.reflection_table['Ih_values'] = flex.double([0.0] * len(self.reflection_table))
    'choose intensities, map to asu, assign unique refl. index'
    self.reflection_table = self.select_optimal_intensities(self.reflection_table)
    self.reflection_table = self.map_indices_to_asu(self.reflection_table)
    #add in a method that automatically updates the weights when the reflection table is changed?
    'assign initial weights (will be statistical weights at this point)'
    self.reflection_table = calc_normE2(self.reflection_table, self.experiments)
    self.reflection_table['wilson_outlier_flag'] = calculate_wilson_outliers(
      self.reflection_table)
    self.weights_for_scaling = self.update_weights_for_scaling(
      self.reflection_table)

  'define a few methods required upon initialisation to set up the data manager'
  def extract_reflections_for_scaling(self, reflection_table, error_model_params=None):
    '''select the reflections with non-zero weight and update scale weights
    object.'''
    n_refl = len(reflection_table)
    weights_for_scaling = self.update_weights_for_scaling(reflection_table)
    sel = weights_for_scaling.get_weights() > 0.0
    #reflections_for_scaling = reflection_table.select(sel)
    sel1 = reflection_table['Esq'] > self.scaling_options['E2min']
    sel2 = reflection_table['Esq'] < self.scaling_options['E2max']
    selection = sel & sel1 & sel2
    reflections_for_scaling = reflection_table.select(selection)
    weights_for_scaling = self.update_weights_for_scaling(reflections_for_scaling,
      error_model_params=error_model_params)
    msg = ('{0} reflections were selected for scale factor determination {sep}'
      'out of {1} reflections. This was based on selection criteria of {sep}'
      'E2min = {2}, E2max = {3}, Isigma_min = {4}, dmin = {5}. {sep}').format(
      reflections_for_scaling.size(), n_refl, self.scaling_options['E2min'],
      self.scaling_options['E2max'], self.scaling_options['Isigma_min'], 
      self.scaling_options['d_min'], sep='\n')
    print(msg)
    return reflections_for_scaling, weights_for_scaling, selection

  def update_weights_for_scaling(self, reflection_table, weights_filter=True, 
                                 error_model_params=None):
    '''set the weights of each reflection to be used in scaling'''
    weights_for_scaling = Weighting(reflection_table)
    print('Updating weights to be used during scaling. \n')
    if weights_filter:
      weights_for_scaling.apply_Isigma_cutoff(reflection_table,
                                              self.scaling_options['Isigma_min'])
      weights_for_scaling.apply_dmin_cutoff(reflection_table,
                                            self.scaling_options['d_min'])
    weights_for_scaling.remove_wilson_outliers(reflection_table)
    if error_model_params:
      weights_for_scaling.apply_aimless_error_model(reflection_table, error_model_params)
    return weights_for_scaling

  def map_indices_to_asu(self, reflection_table):
    '''Create a miller_set object, map to the asu and create a sorted
       reflection table, sorted by asu miller index'''
    u_c = self.experiments.crystal.get_unit_cell().parameters()
    if self.scaling_options['space_group']:
      sg_from_file = self.experiments.crystal.get_space_group().info()
      s_g_symbol = self.scaling_options['space_group']
      crystal_symmetry = crystal.symmetry(unit_cell=u_c, space_group_symbol=s_g_symbol)
      msg = ('WARNING: Manually overriding space group from {0} to {1}. {sep}'
        'If the reflection indexing in these space groups is different, {sep}'
        'bad things may happen!!! {sep}').format(sg_from_file, s_g_symbol, sep='\n')
      print(msg)
    else:
      s_g = self.experiments.crystal.get_space_group()
      crystal_symmetry = crystal.symmetry(unit_cell=u_c, space_group=s_g)
    miller_set = miller.set(crystal_symmetry=crystal_symmetry,
                            indices=reflection_table['miller_index'], anomalous_flag=True)
    miller_set_in_asu = miller_set.map_to_asu()
    reflection_table["asu_miller_index"] = miller_set_in_asu.indices()
    permuted = (miller_set.map_to_asu()).sort_permutation(by_value='packed_indices')
    reflection_table = reflection_table.select(permuted)
    return reflection_table

  def select_optimal_intensities(self, reflection_table):
    '''method to choose which intensities to use for scaling'''
    if (self.scaling_options['integration_method'] == 'sum' or
        self.scaling_options['integration_method'] == 'prf'):
      intstr = self.scaling_options['integration_method']
      reflection_table['intensity'] = (reflection_table['intensity.'+intstr+'.value']
                                       * reflection_table['lp']
                                       / reflection_table['dqe'])
      reflection_table['variance'] = (reflection_table['intensity.'+intstr+'.variance']
                                      * (reflection_table['lp']**2)
                                      / (reflection_table['dqe']**2))
      print(('{0} intensity values will be used for scaling. {sep}').format(
        'Profile fitted' if intstr == 'prf' else 'Summation integrated', sep='\n'))
    #perform a combined prf/sum in a similar fashion to aimless
    elif self.scaling_options['integration_method'] == 'combine':
      int_prf = (reflection_table['intensity.prf.value']
                 * reflection_table['lp'] / reflection_table['dqe'])
      int_sum = (reflection_table['intensity.sum.value']
                 * reflection_table['lp'] / reflection_table['dqe'])
      var_prf = (reflection_table['intensity.prf.variance']
                 * (reflection_table['lp']**2) / (reflection_table['dqe']**2))
      var_sum = (reflection_table['intensity.sum.variance']
                 * (reflection_table['lp']**2) / (reflection_table['dqe']**2))
      Imid = max(int_sum)/2.0
      weight = 1.0/(1.0 + ((int_prf/Imid)**3))
      reflection_table['intensity'] = ((weight * int_prf) + ((1.0 - weight) * int_sum))
      reflection_table['variance'] = ((weight * var_prf) + ((1.0 - weight) * var_sum))
      msg = ('Combined profile/summation intensity values will be used for {sep}'
      'scaling, with an Imid of {0}. {sep}').format(Imid, sep='\n')
      print(msg)
    return reflection_table

  def update_for_minimisation(self, apm):
    '''update the scale factors and Ih for the next iteration of minimisation'''
    basis_fn = self.get_basis_function(apm)
    self.active_derivatives = basis_fn[1]
    self.Ih_table.update_scale_factors(basis_fn[0])
    self.Ih_table.calc_Ih()

  '''define a few methods for saving the data'''
  def save_reflection_table(self, filename):
    ''' Save the reflections to file. '''
    self.reflection_table.as_pickle(filename)

  def save_data_manager(self, filename):
    ''' Save the data manager to file. '''
    data_file = open(filename, 'w')
    pickle.dump(self, data_file)
    data_file.close()


class aimless_Data_Manager(Data_Manager):
  '''Data Manager subclass for implementing XDS parameterisation'''
  def __init__(self, reflections, experiments, scaling_options):
    Data_Manager.__init__(self, reflections, experiments, scaling_options)
    'initialise g-value objects'
    (self.g_absorption, self.g_scale, self.g_decay) = (None, None, None)
    self.g_parameterisation = OrderedDict()
    '''bin reflections, determine outliers, extract reflections and weights for
    scaling and set normalised values.'''
    self.initialise_scale_factors()
    (reflections_for_scaling, weights_for_scaling, selection) = (
      self.extract_reflections_for_scaling(self.reflection_table,
      error_model_params=self.scaling_options['error_model_params']))
    self.Ih_table = single_Ih_table(reflections_for_scaling,
                                    weights_for_scaling.get_weights())
    '''refactor the next two operations into extract_reflections?
    reset the normalised values within the scale_factor object to current'''
    self.g_scale.set_normalised_values(reflections_for_scaling[
      'normalised_rotation_angle'])
    if self.scaling_options['decay_term']:
      self.g_decay.set_normalised_values(reflections_for_scaling[
        'normalised_time_values'])
      self.g_decay.set_d_values(reflections_for_scaling['d'])
    if self.scaling_options['absorption_term']:
      #self.g_absorption.set_values(self.sph_harm_table.select(selection))
      self.g_absorption.set_values(sph_harm_table(reflections_for_scaling,
                                                  self.scaling_options['lmax']))
    print('Completed initialisation of aimless data manager. \n' + '*'*40 + '\n')

  def initialise_scale_factors(self):
    '''initialise scale factors and add to self.active_parameters'''
    print('Initialising scale factor objects. \n')
    self.initialise_scale_term(self.reflection_table)
    self.initialise_decay_term(self.reflection_table)
    self.initialise_absorption_scales(self.reflection_table,
                                      self.scaling_options['lmax'])

  def get_target_function(self, apm):
    '''call the aimless target function method'''
    return target_function(self, apm).return_targets()

  def get_basis_function(self, apm):
    '''call the aimless basis function method'''
    return bf.basis_function(self, apm).return_basis()

  def update_for_minimisation(self, apm):
    '''update the scale factors and Ih for the next iteration of minimisation'''
    basis_fn = self.get_basis_function(apm)
    self.active_derivatives = basis_fn[1]
    self.Ih_table.update_scale_factors(basis_fn[0])
    self.Ih_table.calc_Ih()

  def initialise_decay_term(self, reflection_table):
    '''calculate the relative, normalised rotation angle. Here this is called
    normalised time to allow a different rotation interval compare to the scale
    correction. A SmoothScaleFactor_1D object is then initialised'''
    if self.scaling_options['decay_term']:
      if self.scaling_options['B_factor_interval']:
        rotation_interval = self.scaling_options['B_factor_interval']
      else:
        rotation_interval = 20.0
      rotation_interval = rotation_interval + 0.001
      osc_range = self.experiments.scan.get_oscillation_range()
      if ((osc_range[1] - osc_range[0]) / rotation_interval) % 1 < 0.33:
        #if last bin less than 33% filled'
        n_phi_bins = int((osc_range[1] - osc_range[0]) / rotation_interval)
        'increase rotation interval slightly'
        rotation_interval = (osc_range[1] - osc_range[0])/float(n_phi_bins) + 0.001
      'extend by 0.001 to make sure all datapoints within min/max'
      one_oscillation_width = self.experiments.scan.get_oscillation()[1]
      reflection_table['normalised_time_values'] = ((reflection_table['xyzobs.px.value'].parts()[2]
        * one_oscillation_width) + 0.001)/rotation_interval
      'define the highest and lowest gridpoints: go out two further than the max/min int values'
      #highest_parameter_value = int((max(reflection_table['normalised_time_values'])//1)+3)#was +2
      #lowest_parameter_value = int((min(reflection_table['normalised_time_values'])//1)-2)#was -1
      highest_parameter_value = int((max(reflection_table['normalised_time_values'])//1)+2)#was +2
      lowest_parameter_value = int((min(reflection_table['normalised_time_values'])//1)-1)#was -1
      n_decay_parameters = highest_parameter_value - lowest_parameter_value + 1
      self.g_decay = SF.SmoothScaleFactor_1D_Bfactor(0.0, n_decay_parameters, reflection_table['d'])
      self.g_parameterisation['g_decay'] = self.g_decay
      msg = ('The decay term ScaleFactor object was successfully initialised. {sep}'
        'The B-factor parameter interval has been set to {0} degrees. {sep}'
        '{1} parameters will be used to parameterise the time-dependent decay. {sep}'
        ).format(rotation_interval, n_decay_parameters, sep='\n')
      print(msg)

  def initialise_scale_term(self, reflection_table):
    '''calculate the relative, normalised rotation angle.
    A SmoothScaleFactor_1D object is then initialised'''
    if self.scaling_options['rotation_interval']:
      rotation_interval = self.scaling_options['rotation_interval']
    else:
      rotation_interval = 15.0
    rotation_interval = rotation_interval + 0.001
    osc_range = self.experiments.scan.get_oscillation_range()
    if ((osc_range[1] - osc_range[0])/ rotation_interval) % 1 < 0.33:
      #if last bin less than 33% filled
      n_phi_bins = int((osc_range[1] - osc_range[0])/ rotation_interval)
      'increase rotation interval slightly'
      rotation_interval = (osc_range[1] - osc_range[0])/float(n_phi_bins) + 0.001
    'extend by 0.001 to make sure all datapoints within min/max'
    one_oscillation_width = self.experiments.scan.get_oscillation()[1]
    reflection_table['normalised_rotation_angle'] = ((reflection_table['xyzobs.px.value'].parts()[2]
      * one_oscillation_width) + 0.001)/rotation_interval
    'define the highest and lowest gridpoints: go out two further than the max/min int values'
    #highest_parameter_value = int((max(reflection_table['normalised_rotation_angle'])//1)+3)#was +2
    #lowest_parameter_value = int((min(reflection_table['normalised_rotation_angle'])//1)-2)#was -1
    highest_parameter_value = int((max(reflection_table['normalised_rotation_angle'])//1)+2)
    lowest_parameter_value = int((min(reflection_table['normalised_rotation_angle'])//1)-1)
    n_scale_parameters = highest_parameter_value - lowest_parameter_value + 1
    self.g_scale = SF.SmoothScaleFactor_1D(1.0, n_scale_parameters)
    self.g_scale.set_normalised_values(reflection_table['normalised_rotation_angle'])
    self.g_parameterisation['g_scale'] = self.g_scale
    msg = ('The scale term ScaleFactor object was successfully initialised. {sep}'
      'The scale term parameter interval has been set to {0} degrees. {sep}'
      '{1} parameters will be used to parameterise the time-dependent scale. {sep}'
        ).format(rotation_interval, n_scale_parameters, sep='\n')
    print(msg)

  def initialise_absorption_scales(self, reflection_table, lmax):
    if self.scaling_options['absorption_term']:
      reflection_table = calc_s2d(reflection_table, self.experiments)
      n_abs_params = 0
      for i in range(lmax):
        n_abs_params += (2*(i+1))+1
      self.sph_harm_table = sph_harm_table(reflection_table, lmax)
      self.g_absorption = SF.SphericalAbsorption_ScaleFactor(0.0, n_abs_params,
        self.sph_harm_table)
      self.g_parameterisation['g_absorption'] = self.g_absorption
      msg = ('The absorption term ScaleFactor object was successfully initialised. {sep}'
        'The absorption term will be parameterised by a set of spherical {sep}'
        'harmonics up to an lmax of {0} ({1} parameters). {sep}'
        ).format(lmax, n_abs_params, sep='\n')
      print(msg)
    else:
      reflection_table['phi'] = (reflection_table['xyzobs.px.value'].parts()[2]
                                 * self.experiments.scan.get_oscillation()[1])
    

  def calc_absorption_constraint(self, apm):
    #this should only be called if g_absorption in active params
    idx = apm.active_parameterisation.index('g_absorption')
    start_idx = apm.cumulative_active_params[idx]
    end_idx = apm.cumulative_active_params[idx+1]
    weight = 1e6
    abs_params = apm.x[start_idx:end_idx]
    residual = (weight * (abs_params)**2)
    gradient = (2 * weight * abs_params)
    #need to make sure gradient is returned is same size as gradient calculated in target fn-
    #would only be triggered if refining absorption as same time as another scale factor.
    gradient_vector = flex.double([])
    for i, param in enumerate(apm.active_parameterisation):
      if param != 'g_absorption':
        gradient_vector.extend(flex.double([0.0]*apm.active_params_list[i]))
      elif param == 'g_absorption':
        gradient_vector.extend(gradient)
    return (residual, gradient_vector)

  def normalise_scales_and_B(self):
    if self.scaling_options['decay_term']:
      absorption_scales = self.g_decay.get_scales_of_reflections()
      
      B_values = flex.double(np.log(absorption_scales)) * 2.0 * (self.g_decay.d_values**2)
      #B_parameters = self.g_decay.get_scale_factors()
      B_parameters = self.g_decay.value
      B_new_parameters = B_parameters - flex.double([max(B_values)]*len(B_parameters))
      self.g_decay.update_scale_factors(B_new_parameters)
      #self.g_decay.value = B_new_parameters
      self.g_decay.calculate_scales_and_derivatives()
      absorption_scales = self.g_decay.get_scales_of_reflections()
      B_values = flex.double(np.log(absorption_scales)) * 2.0 * (self.g_decay.d_values**2)
    #scale_factors = self.g_scale.get_scale_factors()
    scales = self.g_scale.get_scales_of_reflections()
    normalised_values = self.g_scale.get_normalised_values()
    first_value = min(normalised_values)
    sel = (normalised_values == first_value)
    initial_scale = list(scales.select(sel))[0]
    #scale_factors = self.g_scale.get_scale_factors()
    scale_factors = self.g_scale.value
    new_scales = scale_factors/initial_scale
    self.g_scale.update_scale_factors(new_scales)
    #self.g_scale.value = new_scales
    self.g_scale.calculate_scales_and_derivatives()

  def expand_scales_to_all_reflections(self):
    expanded_scale_factors = flex.double([1.0]*len(self.reflection_table))
    if not self.scaling_options['multi_mode']:
      self.normalise_scales_and_B()
    "recalculate scales for reflections in sorted_reflection table"
    self.g_scale.set_normalised_values(self.reflection_table['normalised_rotation_angle'])
    self.g_scale.calculate_scales()
    expanded_scale_factors *= self.g_scale.scales
    if self.scaling_options['decay_term']:
      self.g_decay.set_normalised_values(self.reflection_table['normalised_time_values'])
      self.g_decay.set_d_values(self.reflection_table['d'])
      self.g_decay.calculate_scales()
      expanded_scale_factors *= self.g_decay.scales
      absorption_scales = self.g_decay.get_scales_of_reflections()
      B_values = flex.double(np.log(absorption_scales)) * 2.0 * (self.g_decay.d_values**2)
    if self.scaling_options['absorption_term']:
      self.g_absorption.set_values(self.sph_harm_table)
      self.g_absorption.calculate_scales_and_derivatives()
      expanded_scale_factors  *= self.g_absorption.scales
    self.reflection_table['inverse_scale_factor'] = expanded_scale_factors
    print(('Scale factors determined during minimisation have now been applied {sep}'
      'to all reflections. {sep}').format(sep='\n'))
    sel = self.weights_for_scaling.get_weights() != 0.0
    self.reflection_table = self.reflection_table.select(sel)
    #error_model_params = mf.error_scale_LBFGSoptimiser(self.Ih_table, flex.double([1.0,0.123])).x 
    self.weights_for_scaling = self.update_weights_for_scaling(self.reflection_table,
      weights_filter=False, error_model_params=None)
    self.Ih_table = single_Ih_table(self.reflection_table, self.weights_for_scaling.get_weights())
    self.reflection_table['Ih_values'] = self.Ih_table.Ih_table['Ih_values']
    print('A new best estimate for I_h for all reflections has now been calculated. \n')

  def clean_reflection_table(self):
    self.initial_keys.append('inverse_scale_factor')
    self.initial_keys.append('phi')
    for key in self.reflection_table.keys():
      if not key in self.initial_keys:
        del self.reflection_table[key]
    #keep phi column for now for comparing to aimless
    added_columns = ['h_index', 's2', 's2d',
                     'decay_factor', 'angular_scale_factor',
                     'normalised_rotation_angle', 'normalised_time_values',
                     'wilson_outlier_flag', 'centric_flag', 'absorption_factor']
    for key in added_columns:
      del self.reflection_table[key]


class KB_Data_Manager(Data_Manager):
  def __init__(self, reflections, experiments, scaling_options):
    Data_Manager.__init__(self, reflections, experiments, scaling_options)
    self.active_parameters = flex.double([])
    (reflections_for_scaling, weights_for_scaling, selection) = (
      self.extract_reflections_for_scaling(self.reflection_table))
    self.Ih_table = base_Ih_table(reflections_for_scaling, weights_for_scaling.get_weights())
    self.g_parameterisation = OrderedDict()
    if scaling_options['scale_term']:
      self.g_scale = SF.K_ScaleFactor(1.0, n_refl=len(self.Ih_table.Ih_table))
      self.g_parameterisation['g_scale'] = self.g_scale
    if scaling_options['decay_term']:
      self.g_decay = SF.B_ScaleFactor(0.0, reflections_for_scaling['d'])
      self.g_parameterisation['g_decay'] = self.g_decay
    print('Completed initialisation of KB data manager. \n' + '*'*40 + '\n')

  def get_target_function(self, apm):
    '''call the target function method'''
    return target_function_fixedIh(self, apm).return_targets()

  def get_basis_function(self, apm):
    '''call the KB basis function method'''
    return bf.basis_function(self, apm).return_basis()

  def update_for_minimisation(self, apm):
    '''update the scale factors and Ih for the next iteration of minimisation'''
    basis_fn = self.get_basis_function(apm)
    self.active_derivatives = basis_fn[1]
    self.Ih_table.update_scale_factors(basis_fn[0])
    #note - we don't calculate Ih here as using a target instead

  def expand_scales_to_all_reflections(self):
    expanded_scale_factors = flex.double([1.0]*len(self.reflection_table))
    if self.scaling_options['scale_term']:
      self.g_scale.set_n_refl(len(self.reflection_table))
      self.g_scale.calculate_scales()
      expanded_scale_factors *= self.g_scale.scales
    if self.scaling_options['decay_term']:
      self.g_decay.set_d_values(self.reflection_table['d'])
      self.g_decay.calculate_scales()
      expanded_scale_factors *= self.g_decay.scales
    self.reflection_table['inverse_scale_factor'] = expanded_scale_factors
    print(('Scale factors determined during minimisation have now been applied {sep}'
      'to all reflections. {sep}').format(sep='\n'))


class active_parameter_manager(object):
  ''' object to manage the current active parameters during minimisation.
  Separated out to provide a consistent interface between the data manager and
  minimiser.
  Takes in a data manager, needed to access SF objects through g_parameterisation,
  and a param_name list indicating the active parameters.'''
  def __init__(self, Data_Manager, param_name):
    constant_g_values = []
    self.x = flex.double([])
    self.n_active_params = 0
    self.active_parameterisation = []
    self.active_params_list = []
    self.cumulative_active_params = [0]
    for p_type, scalefactor in Data_Manager.g_parameterisation.iteritems():
      if p_type in param_name:
        self.x.extend(scalefactor.get_scale_factors())
        self.n_active_params = len(self.x) #update n_active_params
        self.active_parameterisation.append(p_type)
        n_params = len(scalefactor.get_scale_factors())
        self.active_params_list.append(n_params)
        self.cumulative_active_params.append(self.cumulative_active_params[-1] + n_params)
      else:
        constant_g_values.append(scalefactor.get_scales_of_reflections())
    if len(constant_g_values) > 0.0:
      self.constant_g_values = flex.double(np.prod(np.array(constant_g_values), axis=0))
    else:
      self.constant_g_values = None
    print(('Set up parameter handler for following corrections: {0}\n').format(
      ''.join(i.lstrip('g_')+' ' for i in self.active_parameterisation)))

  def get_current_parameters(self):
    return self.x

class multi_active_parameter_manager(object):
  def __init__(self, Data_Manager, param_name):
    self.apm_list = []
    for DM in Data_Manager.data_managers:
      apm = active_parameter_manager(DM, param_name)
      self.apm_list.append(apm)
    self.active_parameterisation = self.apm_list[0].active_parameterisation
    self.x = flex.double([])
    for apm in self.apm_list:
      self.x.extend(apm.x)
    self.n_active_params_dataset1 = len(self.apm_list[0].x)
    self.n_active_params_dataset2 = len(self.apm_list[1].x)
    self.n_active_params = self.n_active_params_dataset1 + self.n_active_params_dataset2

  def get_current_parameters(self):
    return self.x

class XDS_Data_Manager(Data_Manager):
  '''Data Manager subclass for implementing XDS parameterisation'''
  def __init__(self, reflections, experiments, scaling_options):
    Data_Manager.__init__(self, reflections, experiments, scaling_options)
    #initialise g-value arrays
    (self.g_absorption, self.g_modulation, self.g_decay) = (None, None, None)
    self.g_parameterisation = {}
    self.initialise_scale_factors() #this creates ScaleFactor objects
    (reflections_for_scaling, weights_for_scaling, selection) = (
      self.extract_reflections_for_scaling(self.reflection_table))
    self.Ih_table = single_Ih_table(reflections_for_scaling, weights_for_scaling.get_weights())
    #update normalised values after extracting reflections for scaling
    if scaling_options['modulation']:
      self.g_modulation.set_normalised_values(reflections_for_scaling['normalised_x_values'],
        reflections_for_scaling['normalised_y_values'])
    if scaling_options['decay']:
      self.g_decay.set_normalised_values(reflections_for_scaling['normalised_res_values'],
        reflections_for_scaling['normalised_time_values'])
    if scaling_options['absorption']:
      self.g_absorption.set_normalised_values(reflections_for_scaling['normalised_x_abs_values'],
        reflections_for_scaling['normalised_y_abs_values'],
        reflections_for_scaling['normalised_time_values'])
    print('Completed initialisation of XDS data manager. \n' + '*'*40 + '\n')

  def initialise_scale_factors(self):
    print('Initialising scale factor objects. \n')
    self.bin_reflections_decay()
    if self.scaling_options['absorption']:
      self.bin_reflections_absorption()
    if self.scaling_options['modulation']:
      self.bin_reflections_modulation()

  def get_target_function(self, parameters):
    '''call the xds target function method'''
    if self.scaling_options['parameterization'] == 'log':
      return xds_target_function_log(self, parameters).return_targets()
    else:
      return target_function(self, parameters).return_targets()

  def get_basis_function(self, parameters):
    '''call the xds basis function method'''
    if self.scaling_options['parameterization'] == 'log':
      return bf.xds_basis_function_log(self, parameters).return_basis()
    else:
      return bf.basis_function(self, parameters).return_basis()

  def bin_reflections_decay(self):
    '''bin reflections for decay correction'''
    ndbins = self.scaling_options['n_d_bins']
    if self.scaling_options['rotation_interval']:
      rotation_interval = self.scaling_options['rotation_interval']
    else:
      rotation_interval = 15.0
    osc_range = self.experiments.scan.get_oscillation_range()
    if ((osc_range[1] - osc_range[0]) / rotation_interval) % 1 < 0.33:
      #if last bin less than 33% filled'
      n_phi_bins = int((osc_range[1] - osc_range[0])/ rotation_interval)
      'increase rotation interval slightly'
      rotation_interval = (osc_range[1] - osc_range[0])/float(n_phi_bins) + 0.001
    else:
      rotation_interval = 15.0 + 0.001
      n_phi_bins = int((osc_range[1] - osc_range[0]) / rotation_interval) + 1
    self.n_phi_bins = n_phi_bins
    nzbins = n_phi_bins
    '''Bin the data into resolution and time 'z' bins'''
    zmax = max(self.reflection_table['xyzobs.px.value'].parts()[2]) + 0.001
    zmin = min(self.reflection_table['xyzobs.px.value'].parts()[2]) - 0.001
    resmax = (1.0 / (min(self.reflection_table['d'])**2)) + 0.001
    resmin = (1.0 / (max(self.reflection_table['d'])**2)) - 0.001
    one_resbin_width = (resmax - resmin) / ndbins
    one_time_width = (zmax - zmin) / nzbins
    self.reflection_table['normalised_res_values'] = (
      ((1.0 / (self.reflection_table['d']**2)) - resmin) / one_resbin_width)
    #define the highest and lowest gridpoints: go out one further than the max/min int values
    highest_parameter_value = int((max(self.reflection_table['normalised_res_values'])//1)+2)
    lowest_parameter_value = int((min(self.reflection_table['normalised_res_values'])//1)-1)
    n_res_parameters = highest_parameter_value - lowest_parameter_value + 1
    self.reflection_table['normalised_time_values'] = (
      (self.reflection_table['xyzobs.px.value'].parts()[2] - zmin) / one_time_width)
    if self.scaling_options['decay']:
      #put this here to allow calculation of normalised time values needed for abscor, will rearrange.
      #define the highest and lowest gridpoints: go out one further than the max/min int values
      highest_parameter_value = int((max(self.reflection_table['normalised_time_values'])//1)+2)
      lowest_parameter_value = int((min(self.reflection_table['normalised_time_values'])//1)-1)
      n_time_parameters = highest_parameter_value - lowest_parameter_value + 1
      if self.scaling_options['parameterization'] == 'log':
        self.g_decay = SF.SmoothScaleFactor_2D(0.0, n_res_parameters, n_time_parameters)
      else:
        self.g_decay = SF.SmoothScaleFactor_2D(1.0, n_res_parameters, n_time_parameters)
      self.g_decay.set_normalised_values(self.reflection_table['normalised_res_values'],
        self.reflection_table['normalised_time_values'])
      self.g_parameterisation['g_decay'] = self.g_decay
    msg = ('The decay term ScaleFactor object was successfully initialised. {sep}'
      'The rotational parameter interval has been set to {0} degrees, giving {sep}'
      '{1} phi bins, while the resolution has been binned into {2} bins. {sep}'
      'In total, this ScaleFactor object uses {3} parameters. {sep}'
      ).format(rotation_interval, n_phi_bins, ndbins, 
               n_time_parameters * n_res_parameters, sep='\n')
    print(msg)

  def bin_reflections_modulation(self):
    '''bin reflections for modulation correction'''
    nxbins = nybins = self.scaling_options['n_detector_bins']
    xvalues = self.reflection_table['xyzobs.px.value'].parts()[0]
    (xmax, xmin) = (max(xvalues) + 0.001, min(xvalues) - 0.001)
    yvalues = self.reflection_table['xyzobs.px.value'].parts()[1]
    (ymax, ymin) = (max(yvalues) + 0.001, min(yvalues) - 0.001)
    one_x_width = (xmax - xmin) / float(nxbins)
    one_y_width = (ymax - ymin) / float(nybins)
    self.reflection_table['normalised_x_values'] = ((xvalues - xmin) / one_x_width)
    self.reflection_table['normalised_y_values'] = ((yvalues - ymin) / one_y_width)
    #define the highest and lowest gridpoints: go out one further than the max/min int values
    highest_x_parameter_value = int((max(self.reflection_table['normalised_x_values'])//1)+2)
    lowest_x_parameter_value = int((min(self.reflection_table['normalised_x_values'])//1)-1)
    n_x_parameters = highest_x_parameter_value - lowest_x_parameter_value + 1
    highest_y_parameter_value = int((max(self.reflection_table['normalised_y_values'])//1)+2)
    lowest_y_parameter_value = int((min(self.reflection_table['normalised_y_values'])//1)-1)
    n_y_parameters = highest_y_parameter_value - lowest_y_parameter_value + 1
    if self.scaling_options['parameterization'] == 'log':
      self.g_modulation = SF.SmoothScaleFactor_2D(0.0, n_x_parameters, n_y_parameters)
    else:
      self.g_modulation = SF.SmoothScaleFactor_2D(1.0, n_x_parameters, n_y_parameters)
    self.g_modulation.set_normalised_values(self.reflection_table['normalised_x_values'],
      self.reflection_table['normalised_y_values'])
    self.g_parameterisation['g_modulation'] = self.g_modulation
    msg = ('The modulation term ScaleFactor object was successfully initialised. {sep}'
      'This parameterises a correction across the detector area, using a {sep}'
      'square grid of {0} parameters. {sep}').format(n_x_parameters**2, sep='\n')
    print(msg)


  def bin_reflections_absorption(self):
    '''bin reflections for absorption correction'''
    nxbins = nybins = 4.0
    xvalues = self.reflection_table['xyzobs.px.value'].parts()[0]
    (xmax, xmin) = (max(xvalues) + 0.001, min(xvalues) - 0.001)
    yvalues = self.reflection_table['xyzobs.px.value'].parts()[1]
    (ymax, ymin) = (max(yvalues) + 0.001, min(yvalues) - 0.001)
    one_x_width = (xmax - xmin) / float(nxbins)
    one_y_width = (ymax - ymin) / float(nybins)
    self.reflection_table['normalised_x_abs_values'] = ((xvalues - xmin) / one_x_width)
    self.reflection_table['normalised_y_abs_values'] = ((yvalues - ymin) / one_y_width)
    #define the highest and lowest gridpoints: go out one further than the max/min int values
    highest_x_parameter_value = int((max(self.reflection_table['normalised_x_abs_values'])//1)+1)
    lowest_x_parameter_value = int((min(self.reflection_table['normalised_x_abs_values'])//1))
    n_x_parameters = highest_x_parameter_value - lowest_x_parameter_value + 1
    highest_y_parameter_value = int((max(self.reflection_table['normalised_y_abs_values'])//1)+1)
    lowest_y_parameter_value = int((min(self.reflection_table['normalised_y_abs_values'])//1))
    n_y_parameters = highest_y_parameter_value - lowest_y_parameter_value + 1
    #new code to change time binning to smooth parameters
    highest_parameter_value = int((max(self.reflection_table['normalised_time_values'])//1)+1)
    lowest_parameter_value = int((min(self.reflection_table['normalised_time_values'])//1)-0)
    n_time_parameters = highest_parameter_value - lowest_parameter_value + 1
    #n_time_bins = int((max(self.reflection_table['normalised_time_values'])//1)+1)
    if self.scaling_options['parameterization'] == 'log':
      self.g_absorption = SF.SmoothScaleFactor_GridAbsorption(
        0.0, n_x_parameters, n_y_parameters, n_time_parameters)
    else:
      self.g_absorption = SF.SmoothScaleFactor_GridAbsorption(
        1.0, n_x_parameters, n_y_parameters, n_time_parameters)
    self.g_absorption.set_normalised_values(self.reflection_table['normalised_x_abs_values'],
      self.reflection_table['normalised_y_abs_values'],
      self.reflection_table['normalised_time_values'])
    self.g_parameterisation['g_absorption'] = self.g_absorption
    msg = ('The absorption term ScaleFactor object was successfully initialised. {sep}'
      'This defines a {0} by {0} grid of absorption correction parameters {sep}'
      'at {1} time values, giving a total of {2} parameters. {sep}'
      ).format(n_x_parameters, n_time_parameters, 
               n_time_parameters * (n_x_parameters**2), sep='\n')
    print(msg)


  def bin_reflections_absorption_radially(self):
    from math import pi
    n_rad_bins = 3.0
    n_ang_bins = 8.0
    image_size = self.experiments.detector.to_dict()['panels'][0]['image_size']
    xvalues = self.reflection_table['xyzobs.px.value'].parts()[0]
    yvalues = self.reflection_table['xyzobs.px.value'].parts()[1]
    x_center = image_size[0]/2.0
    y_center = image_size[1]/2.0
    xrelvalues = xvalues - x_center
    yrelvalues = yvalues - y_center
    radial_values = ((xrelvalues**2) + (yrelvalues**2))**0.5
    angular_values = flex.double(np.arctan2(yrelvalues, xrelvalues))
    normalised_radial_values = n_rad_bins * radial_values / (max(radial_values) * 1.001)
    normalised_angular_values = n_ang_bins * (angular_values + pi) / (2.0 * pi)
    #catchers to stop values being exactly on the boundary and causing errors later.
    sel = normalised_angular_values == 8.0
    normalised_angular_values.set_selected(sel, 7.9999)
    sel = normalised_angular_values == 0.0
    normalised_angular_values.set_selected(sel, 0.0001)
    self.reflection_table['normalised_angle_values'] = normalised_angular_values
    self.reflection_table['normalised_radial_values'] = normalised_radial_values
    n_ang_parameters = n_rad_bins + 1
    n_rad_parameters = n_rad_bins + 1
    #new code to change time binning to smooth parameters
    highest_parameter_value = int((max(self.reflection_table['normalised_time_values'])//1)+1)
    lowest_parameter_value = int((min(self.reflection_table['normalised_time_values'])//1)-0)
    n_time_parameters = highest_parameter_value - lowest_parameter_value + 1

  def expand_scales_to_all_reflections(self):
    #currently we have Ih, scales and weights for reflections_for_scaling
    #first calculate the scales for all reflections (Sorted reflections)
    scales_to_expand = []
    if self.scaling_options['modulation']:
      self.g_modulation.set_normalised_values(self.reflection_table['normalised_x_values'],
        self.reflection_table['normalised_y_values'])
      scales_to_expand.append(self.g_modulation.calculate_smooth_scales())
    if self.scaling_options['decay']:
      self.g_decay.set_normalised_values(self.reflection_table['normalised_res_values'],
        self.reflection_table['normalised_time_values'])
      scales_to_expand.append(self.g_decay.calculate_smooth_scales())
    if self.scaling_options['absorption']:
      self.g_absorption.set_normalised_values(self.reflection_table['normalised_x_abs_values'],
        self.reflection_table['normalised_y_abs_values'],
        self.reflection_table['normalised_time_values'])
      scales_to_expand.append(self.g_absorption.calculate_smooth_scales())
    
    if self.scaling_options['parameterization'] == 'log':
      self.reflection_table['inverse_scale_factor'] = flex.double(
        np.exp(np.sum(np.array(scales_to_expand), axis=0)))
    else:
      self.reflection_table['inverse_scale_factor'] = flex.double(
        np.prod(np.array(scales_to_expand), axis=0))
    'remove reflections that were determined as outliers'
    print(('Scale factors determined during minimisation have now been applied {sep}'
      'to all reflections. {sep}').format(sep='\n'))
    self.weights_for_scaling = self.update_weights_for_scaling(
      self.reflection_table, weights_filter=False)
    sel = self.weights_for_scaling.get_weights() != 0.0
    self.reflection_table = self.reflection_table.select(sel)
    self.weights_for_scaling = self.update_weights_for_scaling(
      self.reflection_table, weights_filter=False)
    #now calculate Ih for all reflections.
    self.Ih_table = single_Ih_table(self.reflection_table, 
      self.weights_for_scaling.get_weights())
    print('A new best estimate for I_h for all reflections has now been calculated. \n')

  def clean_reflection_table(self):
    #add keys for additional data that is to be exported
    self.initial_keys.append('inverse_scale_factor')
    for key in self.reflection_table.keys():
      if not key in self.initial_keys:
        del self.reflection_table[key]
    added_columns = ['l_bin_index', 'a_bin_index', 'xy_bin_index', 'h_index',
      'normalised_y_values', 'normalised_x_values', 'normalised_y_abs_values',
      'normalised_x_abs_values', 'normalised_time_values', 
      'normalised_res_values', 'wilson_outlier_flag', 'centric_flag']
    for key in added_columns:
      del self.reflection_table[key]

class multicrystal_datamanager(Data_Manager):
  def __init__(self, reflections1, experiments1, reflections2, experiments2, scaling_options):
    if scaling_options['scaling_method'] == 'xds':
      self.dm1 = XDS_Data_Manager(reflections1, experiments1, scaling_options)
      self.dm2 = XDS_Data_Manager(reflections2, experiments2, scaling_options)
    elif scaling_options['scaling_method'] == 'aimless':
      self.dm1 = aimless_Data_Manager(reflections1, experiments1, scaling_options)
      self.dm2 = aimless_Data_Manager(reflections2, experiments2, scaling_options)
    else:
      assert 0, "Incorrect scaling method passed to multicrystal datamanager (not 'xds' or 'aimless')"
    self.experiments = experiments1 #assume same space group from two json files.?
    self.joined_Ih_table = joined_Ih_table(self.dm1.Ih_table, self.dm2.Ih_table, experiments1)
    self.n_active_params = None
    self.n_active_params_dataset1 = None
    self.n_active_params_dataset2 = None
    self.scaling_options = scaling_options
    self.data_managers = [self.dm1, self.dm2]
    reflections_for_scaling = self.zip_data_together(self.dm1.Ih_table.Ih_table,
                                                     self.dm2.Ih_table.Ih_table)
    weights_for_scaling = reflections_for_scaling['weights']
    self.Ih_table = single_Ih_table(reflections_for_scaling, weights_for_scaling)
    print('Completed initialisation of multicrystal data manager. \n' + '*'*40 + '\n')

  def get_target_function(self, apm):
    '''call the xds target function method'''
    return target_function(self, apm).return_targets()

  def calc_absorption_constraint(self, apm):
    'method only called in aimless scaling'
    R = flex.double([])
    G = flex.double([])
    R.extend(self.dm1.calc_absorption_constraint(apm.apm_list[0])[0])
    R.extend(self.dm2.calc_absorption_constraint(apm.apm_list[1])[0])
    G.extend(self.dm1.calc_absorption_constraint(apm.apm_list[0])[1])
    G.extend(self.dm2.calc_absorption_constraint(apm.apm_list[1])[1])
    return (R, G)

  def zip_data_together(self, Ih_table_1, Ih_table_2):
    joined_reflections = flex.reflection_table()
    h_idx_cumulative_1 = self.joined_Ih_table.h_index_cumulative_array_1
    h_idx_cumulative_2 = self.joined_Ih_table.h_index_cumulative_array_2
    for i in range(len(h_idx_cumulative_1)-1):
      joined_reflections.extend(Ih_table_1[h_idx_cumulative_1[i]:
                                           h_idx_cumulative_1[i+1]])
      joined_reflections.extend(Ih_table_2[h_idx_cumulative_2[i]:
                                           h_idx_cumulative_2[i+1]])
    return joined_reflections

  def zip_together_scales(self, scales1, scales2):
    scales_1_expanded = scales1 * self.joined_Ih_table.h_expand_mat_1
    scales_2_expanded = scales2 * self.joined_Ih_table.h_expand_mat_2
    return scales_1_expanded + scales_2_expanded

  def zip_together_derivatives(self, apm, derivs1, derivs2):
    n_refl_1 = len(self.dm1.Ih_table.Ih_table)
    n_refl_2 = len(self.dm2.Ih_table.Ih_table)
    n_param_1 = apm.n_active_params_dataset1
    n_param_2 = apm.n_active_params_dataset2
    active_derivatives = sparse.matrix((n_refl_1 + n_refl_2),(n_param_1 + n_param_2))                                   
    derivs_1_T = derivs1.transpose()
    expanded_1 = derivs_1_T * self.joined_Ih_table.h_expand_mat_1
    active_derivatives.assign_block(expanded_1.transpose(), 0, 0)
    derivs_2_T = derivs2.transpose()
    expanded_2 = derivs_2_T * self.joined_Ih_table.h_expand_mat_2
    active_derivatives.assign_block(expanded_2.transpose(), 0, n_param_1)
    return active_derivatives

  def update_for_minimisation(self, apm):
    '''update the scale factors and Ih for the next iteration of minimisation,
    update the x values from the amp to the individual apms, as this is where 
    basis functions, target functions etc get access to the parameters.'''
    apm.apm_list[0].x = apm.x[:apm.n_active_params_dataset1]
    apm.apm_list[1].x = apm.x[apm.n_active_params_dataset1:]
    basis_fn_1 = self.dm1.get_basis_function(apm.apm_list[0])
    basis_fn_2 = self.dm2.get_basis_function(apm.apm_list[1])
    self.dm1.Ih_table.Ih_table['inverse_scale_factor'] = basis_fn_1[0]
    self.dm2.Ih_table.Ih_table['inverse_scale_factor'] = basis_fn_2[0]
    self.active_derivatives = self.zip_together_derivatives(apm, basis_fn_1[1], basis_fn_2[1])
    self.Ih_table.Ih_table['inverse_scale_factor'] = self.zip_together_scales(
      basis_fn_1[0], basis_fn_2[0])
    self.Ih_table.calc_Ih()

  def expand_scales_to_all_reflections(self):
    self.dm1.expand_scales_to_all_reflections()
    self.dm2.expand_scales_to_all_reflections()

  def join_multiple_datasets(self):
    self.joined_Ih_table = joined_Ih_table(self.dm1.Ih_table, self.dm2.Ih_table, self.experiments)
    joined_reflections = flex.reflection_table()
    h_idx_cumulative_1 = self.joined_Ih_table.h_index_cumulative_array_1
    h_idx_cumulative_2 = self.joined_Ih_table.h_index_cumulative_array_2
    for i in range(len(h_idx_cumulative_1)-1):
      joined_reflections.extend(self.dm1.reflection_table[h_idx_cumulative_1[i]:
                                                          h_idx_cumulative_1[i+1]])
      joined_reflections.extend(self.dm2.reflection_table[h_idx_cumulative_2[i]:
                                                          h_idx_cumulative_2[i+1]])
    self.reflection_table = joined_reflections
    self.weights_for_scaling = Weighting(self.reflection_table)
    #self.reflection_table['wilson_outlier_flag'] = (
    #calculate_wilson_outliers(self.reflection_table, self.experiments))
    #self.weights_for_scaling.remove_wilson_outliers(self.reflection_table)
    self.Ih_table = single_Ih_table(self.reflection_table, self.weights_for_scaling.get_weights())
    self.reflection_table['Ih_values'] = self.Ih_table.Ih_table['Ih_values']

class targeted_datamanager(Data_Manager):
  def __init__(self, reflections1, experiments1, reflections_scaled, scaling_options):
    #first assume that the Ih_values of reflections_scaled are the best estimates
    osc_range = experiments1.scan.get_oscillation_range()
    self.experiments = experiments1
    self.scaling_options = scaling_options
    if osc_range[1]-osc_range[0]<1000.0: 
      #usually would have it osc_range <10, big here just for testing on LCY dataset
      'first make a simple KB data manager'
      self.dm1 = KB_Data_Manager(reflections1, experiments1, scaling_options)
      'now extract Ih values from reflections_scaled'
      self.target_dm = Data_Manager(reflections_scaled, experiments1, scaling_options)
      (target_reflections, target_weights, target_selection) = (
      self.extract_reflections_for_scaling(self.target_dm.reflection_table))
      target_Ih_table = single_Ih_table(target_reflections, target_weights.get_weights())
      'find common reflections in the two datasets'
      for i, miller_idx in enumerate(self.dm1.Ih_table.Ih_table['asu_miller_index']):
        sel = target_Ih_table.Ih_table['asu_miller_index'] == miller_idx
        Ih_values = target_Ih_table.Ih_table['Ih_values'].select(sel)
        if Ih_values:
          self.dm1.Ih_table.Ih_table['Ih_values'][i] = Ih_values[0]
          #select [0] above as all Ih values are the same for asu_miller_idx
        else:
          self.dm1.Ih_table.Ih_table['Ih_values'][i] = 0.0
          #should already be zero but set again just in case
      'select only those reflections matched in the target dataset'
      sel = self.dm1.Ih_table.Ih_table['Ih_values'] != 0.0
      self.dm1.Ih_table.Ih_table = self.dm1.Ih_table.Ih_table.select(sel)
      'update the data in the SF objects'
      if self.scaling_options['decay_term']:
        self.dm1.g_decay.set_d_values(self.dm1.g_decay.d_values.select(sel))
      if self.scaling_options['scale_term']:
        self.dm1.g_scale.set_n_refl(len(self.dm1.Ih_table.Ih_table))
      self.g_parameterisation = self.dm1.g_parameterisation
    else:
      #do full aimless/xds scaling of one dataset against the other?
      #self.target_refl_table = reflections_scaled
      assert 0, "no methods specified yet for scaling a large dataset against another"

  def get_target_function(self, apm):
    '''call the target function method'''
    return self.dm1.get_target_function(apm)

  def get_basis_function(self, apm):
    '''call the KB basis function method'''
    return self.dm1.get_basis_function(apm)

  def update_for_minimisation(self, apm):
    '''update the scale factors and Ih for the next iteration of minimisation'''
    return self.dm1.update_for_minimisation(apm)

  def expand_scales_to_all_reflections(self):
    return self.dm1.expand_scales_to_all_reflections()
