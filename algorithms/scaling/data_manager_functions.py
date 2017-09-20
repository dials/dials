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
from data_quality_assessment import R_meas, R_pim

class Data_Manager(object):
  '''Data Manager takes a params parsestring containing the parsed
     integrated.pickle and integrated_experiments.json files'''
  def __init__(self, reflections, experiments, scaling_options):
    'General attributes relevant for all parameterisations'
    self.experiments = experiments
    self.reflection_table = reflections[0]
    self.initial_keys = [key for key in self.reflection_table.keys()]
    self.reflection_table['x_value'] = self.reflection_table['xyzobs.px.value'].parts()[0]
    self.reflection_table['y_value'] = self.reflection_table['xyzobs.px.value'].parts()[1]
    self.reflection_table['z_value'] = self.reflection_table['xyzobs.px.value'].parts()[2]
    self.reflection_table['inverse_scale_factor'] = flex.double([1.0]*len(self.reflection_table))
    self.reflection_table['Ih_values'] = flex.double([0.0]*len(self.reflection_table))
    self.sorted_by_miller_index = False
    self.sorted_reflections = None
    self.filtered_reflections = None
    self.Ih_array = None
    self.h_index_counter_array = None
    self.h_index_cumulative_array = None
    self.scaling_options = scaling_options
    #set choice of intensity and filter out bad values of d, variance etc.
    #fltering methods could be replaced by use of reflection table flags?
    self.select_optimal_intensities(scaling_options)
    #sort the reflection table
    self.filter_and_sort_reflections()
    self.weights_for_scaling = Weighting(self.sorted_reflections)
    self.weights_for_scaling.apply_Isigma_cutoff(self.sorted_reflections,
                                                 scaling_options['Isigma_min'])
    self.weights_for_scaling.apply_dmin_cutoff(self.sorted_reflections, 
                                               scaling_options['d_min'])

  def filter_and_sort_reflections(self):
    self.filtered_reflections = copy.deepcopy(self.reflection_table)
    self.filter_negative_variances()
    self.filter_data('d', -1.0, 0.0)
    self.map_indices_to_asu()
    self.assign_h_index()

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

  def select_optimal_intensities(self, scaling_options):
    if (scaling_options['integration_method'] == 'sum' or
        scaling_options['integration_method'] == 'prf'):
      intstr = scaling_options['integration_method']
      self.reflection_table['intensity'] = (self.reflection_table['intensity.'+intstr+'.value']
                                            * self.reflection_table['lp']
                                            / self.reflection_table['dqe'])
      self.reflection_table['variance'] = (self.reflection_table['intensity.'+intstr+'.variance']
                                           * (self.reflection_table['lp']**2)
                                           / (self.reflection_table['dqe']**2))
    elif scaling_options['integration_method'] == 'combine':
      scaling_options['integration_method'] = 'prf'
      self.select_optimal_intensities(scaling_options)
      '''self.reflection_table['intensity'] = self.reflection_table['intensity.prf.value']
      self.reflection_table['variance'] = self.reflection_table['intensity.prf.variance']
      self.filter_and_sort_reflections()
      Iraw_mean = flex.mean(self.sorted_reflections['intensity.sum.value'])
      for Imid in [0.4 * Iraw_mean, 0.6 * Iraw_mean, 0.8 * Iraw_mean, 1.1 * Iraw_mean, 1.2 * Iraw_mean]:
        weights = 1.0/(1.0 + ((self.sorted_reflections['intensity.sum.value']/Imid)**3.0))
        I = ((weights * self.sorted_reflections['intensity.prf.value'])
            + ((1.0 - weights) * self.sorted_reflections['intensity.sum.value']))
        self.sorted_reflections['intensity'] = I
        sigma = ((weights * (self.sorted_reflections['intensity.prf.variance']**0.5))
           + ((1.0 - weights) * (self.sorted_reflections['intensity.sum.variance']**0.5)))
        self.calc_Ih()
        print R_meas(self)'''

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
        if (self.sorted_reflections['asu_miller_index'][i] ==
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

  '''define a number of data filtering options'''
  def filter_data(self, reflection_table_key, lower, upper):
    '''Filter reflection data for a given measurement variable and limits'''
    bad_data = select_variables_in_range(
      self.filtered_reflections[reflection_table_key], lower, upper)
    inv_sel = ~bad_data
    self.filtered_reflections = self.filtered_reflections.select(inv_sel)

  def filter_I_sigma(self, ratio):
    sel = flex.bool()
    for index, intensity in enumerate(self.filtered_reflections['intensity']):
      if intensity/(self.filtered_reflections['variance'][index]**0.5) > ratio:
        sel.append(True)
      else:
        sel.append(False)
    self.filtered_reflections = self.filtered_reflections.select(sel)

  def filter_negative_intensities(self):
    '''return boolean selection of a given variable range'''
    sel = flex.bool()
    for intensity in self.filtered_reflections['intensity']:
      if intensity <= 0.0:
        sel.append(False)
      else:
        sel.append(True)
    self.filtered_reflections = self.filtered_reflections.select(sel)

  def filter_negative_variances(self):
    '''return boolean selection of a given variable range'''
    sel = flex.bool()
    for intensity in self.filtered_reflections['variance']:
      if intensity <= 0.0:
        sel.append(False)
      else:
        sel.append(True)
    self.filtered_reflections = self.filtered_reflections.select(sel)


  '''define a few methods for saving the data'''
  def save_sorted_reflections(self, filename):
    ''' Save the reflections to file. '''
    self.sorted_reflections.as_pickle(filename)

  def save_data_manager(self, filename):
    ''' Save the data manager to file. '''
    data_file = open(filename, 'w')
    pickle.dump(self, data_file)
    data_file.close()

  def calc_Ih(self):
    '''calculate the current best estimate for I for each reflection group'''
    intensities = self.sorted_reflections['intensity']
    variances = self.sorted_reflections['variance']
    scale_factors = self.sorted_reflections['inverse_scale_factor']
    gsq = (((scale_factors)**2) / variances)
    sumgsq = flex.double(np.add.reduceat(gsq, self.h_index_cumulative_array[:-1]))
    gI = ((scale_factors * intensities) / variances)
    sumgI = flex.double(np.add.reduceat(gI, self.h_index_cumulative_array[:-1]))
    self.Ih_array = sumgI / sumgsq
    self.sorted_reflections['Ih_values'] = flex.double(
      np.repeat(self.Ih_array, self.h_index_counter_array))


class Weighting(object):
  def __init__(self, reflection_table):
    '''set initial weighting to be a statistical weighting'''
    self.scale_weighting = 1.0/reflection_table['variance']

  def get_weights(self):
    return self.scale_weighting

  def set_unity_weighting(self, reflection_table):
    self.scale_weighting = flex.double([1.0]*len(reflection_table['variance']))

  def apply_Isigma_cutoff(self, reflection_table, ratio):
    sel = flex.bool()
    for i, intensity in enumerate(reflection_table['intensity']):
      if ratio > intensity/(reflection_table['variance'][i]**0.5):
        sel.append(True)
      else:
        sel.append(False)
    self.scale_weighting.set_selected(sel, 0.0)

  def apply_dmin_cutoff(self, reflection_table, d_cutoff_value):
    sel = flex.bool()
    for i, d in enumerate(reflection_table['d']):
      if d_cutoff_value > d:
        sel.append(True)
      else:
        sel.append(False)
    self.scale_weighting.set_selected(sel, 0.0)
    print len(self.scale_weighting) - self.scale_weighting.count(0.0)


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
    self.active_bin_index = None
    self.active_param_size = None
    self.constant_g_values = None
    self.initialise_scale_factors()

  def initialise_scale_factors(self):
    self.bin_reflections_decay()
    self.bin_reflections_absorption_radially()
    self.bin_reflections_modulation()
    param_type = 'log'
    if param_type == 'log':
      self.g_parameterisation = {
        'g_absorption' : {'index': 'a_bin_index', 'parameterisation' : self.g_absorption},
        'g_modulation' : {'index': 'xy_bin_index', 'parameterisation' : self.g_modulation},
        'g_decay' : {'index': 'l_bin_index', 'parameterisation' : self.g_decay}}
    else:
      self.g_parameterisation = {
        'g_absorption' : {'index': 'a_bin_index', 'parameterisation' : self.g_absorption},
        'g_modulation' : {'index': 'xy_bin_index', 'parameterisation' : self.g_modulation},
        'g_decay' : {'index': 'l_bin_index', 'parameterisation' : self.g_decay}}

  def get_target_function(self):
    '''call the xds target function method'''
    return xds_target_function(self).return_targets()

  def get_basis_function(self, parameters):
    '''call the xds basis function method'''
    return xds_basis_function(self, parameters).return_basis()

  def update_for_minimisation(self, parameters):
    '''update the scale factors and Ih for the next iteration of minimisation'''
    self.sorted_reflections['inverse_scale_factor'] = self.get_basis_function(parameters)[0]
    self.calc_Ih()

  def bin_reflections_decay(self):
    '''bin reflections for decay correction'''
    ndbins = self.binning_parameters['n_d_bins']
    nzbins = self.binning_parameters['n_z_bins']
    '''Bin the data into resolution and time 'z' bins'''
    zmax = max(self.filtered_reflections['z_value'])
    zmin = min(self.filtered_reflections['z_value'])
    resmax = (1.0 / (min(self.filtered_reflections['d'])**2))
    resmin = (1.0 / (max(self.filtered_reflections['d'])**2))
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
    if firstbin_index.count(-1) > 0 or secondbin_index.count(-1) > 0:
      raise ValueError('Unable to bin data for decay in scaling initialisation')
    self.sorted_reflections['l_bin_index'] = (firstbin_index
                          + (secondbin_index * ndbins))
    self.bin_boundaries = {'d' : d_bins, 'z_value' : z_bins}
    self.g_decay = flex.double([1.0] * (ndbins * nzbins))

  def bin_reflections_modulation(self):
    '''bin reflections for modulation correction'''
    nxbins = nybins = self.binning_parameters['n_detector_bins']
    xvalues = self.sorted_reflections['x_value']
    (xmax, xmin) = (max(xvalues), min(xvalues))
    yvalues = self.sorted_reflections['y_value']
    (ymax, ymin) = (max(yvalues), min(yvalues))
    x_bins = (((flex.double(range(0, nxbins + 1)) * (xmax - xmin) / (nxbins)))
              + flex.double([xmin] * (nxbins + 1)))
    x_bins[0] = x_bins[0] - 0.001
    y_bins = ((flex.double(range(0, nybins + 1)) * (ymax - ymin) /(nybins))
              + flex.double([ymin] * (nybins + 1)))
    y_bins[0] = y_bins[0] - 0.001
    firstbin_index = flex.int([-1] * len(self.sorted_reflections['z_value']))
    secondbin_index = flex.int([-1] * len(self.sorted_reflections['z_value']))
    for i in range(nxbins):
      selection = select_variables_in_range(self.sorted_reflections['x_value'],
                                            x_bins[i], x_bins[i+1])
      firstbin_index.set_selected(selection, i)
    for i in range(nybins):
      selection = select_variables_in_range(self.sorted_reflections['y_value'],
                                            y_bins[i], y_bins[i+1])
      secondbin_index.set_selected(selection, i)
    if firstbin_index.count(-1) > 0 or secondbin_index.count(-1) > 0:
      raise ValueError('Unable to bin data for modulation in scaling initialisation')
    self.sorted_reflections['xy_bin_index'] = (firstbin_index
                                               + (secondbin_index * nxbins))
    self.g_modulation = flex.double([1.0] * (nxbins ** 2))

  def bin_reflections_absorption(self):
    '''bin reflections for absorption correction'''
    nxbins = nybins = self.binning_parameters['n_absorption_positions']
    nzbins = self.binning_parameters['n_z_bins']
    '''Bin the data into detector position and time 'z' bins'''
    z_bins = self.bin_boundaries['z_value']
    #define simple detector area map#
    xvalues = self.sorted_reflections['x_value']
    (xmax, xmin) = (max(xvalues), min(xvalues))
    yvalues = self.sorted_reflections['y_value']
    (ymax, ymin) = (max(yvalues), min(yvalues))
    x_bins = ((flex.double(range(0, nxbins + 1)) * (xmax - xmin) / (nxbins))
              + flex.double([xmin] * (nxbins + 1)))
    x_bins[0] = x_bins[0]-0.001
    y_bins = ((flex.double(range(0, nybins + 1)) * (ymax - ymin) / (nybins))
              + flex.double([ymin] * (nybins + 1)))
    y_bins[0] = y_bins[0]-0.001
    firstbin_index = flex.int([-1] * len(self.sorted_reflections['z_value']))
    secondbin_index = flex.int([-1] * len(self.sorted_reflections['z_value']))
    for i in range(nxbins):
      selection1 = select_variables_in_range(self.sorted_reflections['x_value'],
                                             x_bins[i], x_bins[i+1])
      for j in range(nybins):
        selection2 = select_variables_in_range(self.sorted_reflections['y_value'],
                                               y_bins[j], y_bins[j+1])
        firstbin_index.set_selected(selection1 & selection2,
                                    ((i * nybins) + j))
    for i in range(nzbins):
      selection = select_variables_in_range(self.sorted_reflections['z_value'],
                                            z_bins[i], z_bins[i+1])
      secondbin_index.set_selected(selection, i)
    if firstbin_index.count(-1) > 0 or secondbin_index.count(-1) > 0:
      raise ValueError('Unable to bin data for absorption in scaling initialisation')
    self.sorted_reflections['a_bin_index'] = (firstbin_index
                                              + (secondbin_index * nxbins * nybins))
    self.g_absorption = flex.double([1.0] * ((nxbins**2) * nzbins))

  def bin_reflections_absorption_radially(self):
    '''bin reflections for absorption correction'''
    from math import pi
    nxbins = nybins = self.binning_parameters['n_absorption_positions']
    nzbins = self.binning_parameters['n_z_bins']
    '''Bin the data into detector position and time 'z' bins'''
    z_bins = self.bin_boundaries['z_value']
    #define simple detector area map#
    xvalues = self.sorted_reflections['x_value']
    (xmax, xmin) = (max(xvalues), min(xvalues))
    yvalues = self.sorted_reflections['y_value']
    (ymax, ymin) = (max(yvalues), min(yvalues))
    xrelvalues = xvalues - ((xmax - xmin) / 2.0) #!may need better definition of centerpoint
    yrelvalues = yvalues - ((ymax - ymin) / 2.0) #!may need better definition of centerpoint
    radial_bins = [0.0, ymax / 6.0, 2.0 * ymax / 6.0, (2.0**0.5) * ymax / 2.0]
    angular_bins = [0, pi / 4.0, 2.0 * pi / 4.0, 3.0 * pi / 4.0, pi,
                    5.0 * pi / 4.0, 6.0 * pi / 4.0, 7.0 * pi / 4.0, 2.0 * pi]
    radial_values = ((xrelvalues**2) + (yrelvalues**2))**0.5
    angular_values = np.arccos(yrelvalues/radial_values)
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
    self.sorted_reflections['a_bin_index'] = (firstbin_index + (secondbin_index
                                              * (len(angular_bins) - 1)
                                              * (len(radial_bins) - 1)))
    self.g_absorption = flex.double([1.0] * ((len(angular_bins)-1)
                                             * (len(radial_bins)-1) * nzbins))

  def bin_reflections_absorption_smartly(self):
    import gridding_entropy_dm as gedm
    min_in_each_abs_bin = 30

    nzbins = self.binning_parameters['n_z_bins']
    '''Bin the data into detector position and time 'z' bins'''
    z_bins = self.bin_boundaries['z_value']
    #define simple detector area map#
    xvalues = self.sorted_reflections['x_value']
    (xmax, xmin) = (max(xvalues), min(xvalues))
    yvalues = self.sorted_reflections['y_value']
    (ymax, ymin) = (max(yvalues), min(yvalues))

    xrelvalues = self.sorted_reflections['x_value']*6.0/(xmax+0.001)
    yrelvalues = self.sorted_reflections['y_value']*6.0/(ymax+0.001)

    secondbin_index = flex.int([-1] * len(self.sorted_reflections['z_value']))
    for i in range(nzbins):
      selection = select_variables_in_range(
        self.sorted_reflections['z_value'], z_bins[i], z_bins[i+1])
      secondbin_index.set_selected(selection, i)

    bin_boundaries_x1 = [0, 2, 4, 6]
    bin_boundaries_y1 = [0, 2, 4, 6]
    reflections_for_gridding = gedm.refl_table(xrelvalues, yrelvalues, 
      secondbin_index, nzbins, self.h_index_counter_array, self.h_index_cumulative_array )
    print gedm.calc_entropy_2D(reflections_for_gridding, bin_boundaries_x1, bin_boundaries_y1)
    optimal_boundaries = gedm.perform_optimal_divide(reflections_for_gridding, 
      bin_boundaries_x1, bin_boundaries_y1, min_in_each_abs_bin)
    print len(optimal_boundaries)
    
    firstbin_index = flex.int([-1] * len(self.sorted_reflections['z_value']))
    for i in range(len(optimal_boundaries)):
      selection1 = select_variables_in_range(xrelvalues, optimal_boundaries[i][0][0], optimal_boundaries[i][0][1])
      selection2 = select_variables_in_range(yrelvalues, optimal_boundaries[i][1][0], optimal_boundaries[i][1][1])
      firstbin_index.set_selected(selection1 & selection2, i)
    
    if firstbin_index.count(-1) > 0 or secondbin_index.count(-1) > 0:
      raise ValueError('Unable to bin data for absorption in scaling initialisation')
    self.sorted_reflections['a_bin_index'] = (firstbin_index
                          + (secondbin_index * len(optimal_boundaries)))
    self.g_absorption = flex.double([1.0] * (len(optimal_boundaries) * nzbins))


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
    added_columns = ['l_bin_index','a_bin_index', 'xy_bin_index', 'h_index', 
                      'asu_miller_index']
    for key in added_columns:
      del self.sorted_reflections[key]  

  def reject_outliers(self, tolerance, niter):
    '''Identify outliers using the method of aimless - needs some work'''
    for _ in range(niter):
      Ihl = self.sorted_reflections['intensity']
      variances = self.sorted_reflections['variance']
      scale_factors = self.sorted_reflections['inverse_scale_factor']
      Good_reflections = flex.bool([True]*len(self.sorted_reflections['d']))
      print len(Good_reflections)
      for h in range(len(self.h_index_counter_array)):     
        lsum = self.h_index_counter_array[h]
        if lsum == 2:
          indexer = self.h_index_cumulative_array[h]
          delta1 = (Ihl[indexer] - (scale_factors[indexer] * Ihl[indexer + 1])
              / ((variances[indexer] + ((scale_factors[indexer]**2)
                            * variances[indexer + 1]))**0.5))
          delta2 = (Ihl[indexer + 1] - (scale_factors[indexer + 1] * Ihl[indexer])
              / ((variances[indexer + 1] + ((scale_factors[indexer + 1]**2)
                              * variances[indexer]))**0.5))
          if (abs(delta1/Ihl[indexer + 1]) > tolerance or 
              abs(delta2/Ihl[indexer])  > tolerance):
            #if delta1 > 5.0:
            #    print (delta1, Ihl[indexer + 1])
            #if delta2 > 5.0:
            #    print (delta2, Ihl[indexer])
            Good_reflections[indexer] = False
            Good_reflections[indexer+1] = False
            #outlier_list.append(indexer)
            #outlier_list.append(indexer + 1) #reject both if large difference
        if lsum > 2:
          delta_list=[]
          for i in range(lsum):
            sumgI = 0.0
            sumgsq = 0.0
            I_others_list = np.array([])
            for j in range(lsum):
              if j != i:
                indexer = j + self.h_index_cumulative_array[h]
                sumgI += (Ihl[indexer]*scale_factors[indexer] / variances[indexer])
                sumgsq += (scale_factors[indexer]**2 / variances[indexer])
                I_others_list = np.append(I_others_list,Ihl[indexer])
            Ih_others = sumgI / sumgsq
            indexer = i + self.h_index_cumulative_array[h]
            #print list(I_others_list)
            others_std = np.std(I_others_list, ddof=1)
            delta = ((Ihl[indexer] - (scale_factors[indexer] * Ih_others))
                / ((variances[indexer] + ((scale_factors * others_std)**2))**0.5))
            delta_list.append(delta) 
          max_delta = max([abs(x[0]) for x in delta_list])
          if max_delta/Ih_others > tolerance:
            #if max_delta  > 5.0:
            #    print (max_delta, Ih_others)
            ngreater = len([i for i in delta_list if i > 0])
            nless = len([i for i in delta_list if i < 0])
            if ngreater == 1:
              for n, m in enumerate(delta_list):
                if m > 0:
                  Good_reflections[n + self.h_index_cumulative_array[h]] = False
            if nless == 1:
              for n, m in enumerate(delta_list):
                if m < 0:
                  Good_reflections[n + self.h_index_cumulative_array[h]] = False
            else:
              n = delta_list.index(max_delta)
              Good_reflections[n + self.h_index_cumulative_array[h]] = False
      #print outlier_list  
      print "number of outliers = "+ str(Good_reflections.count(False))
      self.sorted_reflections = self.sorted_reflections.select(Good_reflections)
      self.assign_h_index()
      g1values = flex.double([self.g_decay[i]
        for i in self.sorted_reflections['l_bin_index']])
      g2values = flex.double([self.g_absorption[i]
        for i in self.sorted_reflections['a_bin_index']])
      g3values = flex.double([self.g_modulation[i]
        for i in self.sorted_reflections['xy_bin_index']])
      self.sorted_reflections['inverse_scale_factor'] = (g1values 
        * g2values * g3values)
      if Good_reflections.count(False) == 0:
        break 
    #return Good_reflections


def select_variables_in_range(variable_array, lower_limit, upper_limit):
  '''return boolean selection of a given variable range'''
  sel = flex.bool()
  for variable in variable_array:
    if lower_limit < variable <= upper_limit:
      sel.append(True)
    else:
      sel.append(False)
  return sel
