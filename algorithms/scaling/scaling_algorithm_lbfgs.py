from dials.array_family import flex
from dials.util.options import OptionParser
import minimiser_functions as mf
import data_manager_functions as dmf
from data_plotter import (plot_data, save_data, plot_data_absorption, 
                          plot_data_modulation)
from data_quality_assessment import R_meas, R_pim
import matplotlib.pyplot as plt
import time

def scaling_lbfgs(inputparams, ndbins, nzbins, n_absorption_positions, 
                  n_detector_bins, niter=2,integration_method='summation_integration',
                  decay_correction_rescaling=False,
                  correction_order=['absorption','modulation','decay']):
    """This algorithm loads a reflection table and json file and
    performs a Kabsch-scaling parameterisation using an lbfgs minimiser
    """
    if integration_method == 'profile_fitting':
        int_str = str('intensity.prf.value')
        var_str = str('intensity.prf.variance')
    elif integration_method == 'summation_integration':
        int_str = str('intensity.sum.value')
        var_str = str('intensity.sum.variance')
    else:
        raise ValueError('Incorrect choice of integration_method')

    loaded_reflections = dmf.Data_Manager(inputparams,int_str,var_str)

    """Filter out zero/negative values of d, variance, etc"""
    #loaded_reflections.filter_data(int_str, -10.0, 0.0)
    loaded_reflections.filter_data(var_str, -10.0, 0.0)
    loaded_reflections.filter_data('d', -1.0, 0.0)

    '''Map the indices to the asu and also sort the reflection table by miller index'''
    loaded_reflections.map_indices_to_asu()

    '''determine gridding index for scale parameters '''
    loaded_reflections.bin_reflections_dz(ndbins=ndbins, nzbins=nzbins)
    loaded_reflections.bin_reflections_absorption(npos=n_absorption_positions)
    loaded_reflections.bin_reflections_modulation(ngridpoints=n_detector_bins)
    
    '''assign a unique reflection index to each group of reflections'''
    loaded_reflections.assign_h_index()
    loaded_reflections.create_index_tuple()
    loaded_reflections.scale_by_LP_and_dqe()

    '''call the optimiser on the Data Manager object''' 
    minimised = mf.optimiser(loaded_reflections, correction=correction_order[0],
                             decay_correction_rescaling=decay_correction_rescaling)
    minimised = mf.optimiser(minimised.data_manager, correction=correction_order[1],
                             decay_correction_rescaling=decay_correction_rescaling)
    minimised = mf.optimiser(minimised.data_manager, correction=correction_order[2],
                             decay_correction_rescaling=decay_correction_rescaling)
    
    '''repeat n times - replace with convergence criteria and max iter instead?'''
    for i in range(0,niter):
        minimised = mf.optimiser(minimised.data_manager, correction=correction_order[0],
                             decay_correction_rescaling=decay_correction_rescaling)
        minimised = mf.optimiser(minimised.data_manager, correction=correction_order[1],
                             decay_correction_rescaling=decay_correction_rescaling)
        minimised = mf.optimiser(minimised.data_manager, correction=correction_order[2],
                             decay_correction_rescaling=decay_correction_rescaling)
    'fill in gvalues, Ih_values to reflection table'
    minimised.data_manager.sorted_reflections['Ih'] = flex.double([minimised.data_manager.Ih_array[i] for i in
        minimised.data_manager.sorted_reflections['h_index']])
    minimised.data_manager.sorted_reflections['g_decay'] = flex.double([minimised.data_manager.g_values[i] for i in
        minimised.data_manager.sorted_reflections['l_bin_index']])
    minimised.data_manager.sorted_reflections['g_absorption'] = flex.double([minimised.data_manager.g2_values[i] for i in
        minimised.data_manager.sorted_reflections['a_bin_index']])
    minimised.data_manager.sorted_reflections['g_modulation'] = flex.double([minimised.data_manager.g3_values[i] for i in
        minimised.data_manager.sorted_reflections['xy_bin_index']])
    
    return minimised


if __name__ == "__main__":
    '''load datafiles - can be replaced with sys.argv commands'''
    filepath = "/Users/whi10850/Documents/test_data/integrate/13_integrated.pickle"
    json_filepath = "/Users/whi10850/Documents/test_data/integrate/13_integrated_experiments.json"
    parsestring = OptionParser(read_experiments=True, read_reflections=True,
                               check_format=False)
    params, options = parsestring.parse_args([json_filepath, filepath])

    '''do minimisation of g-factors'''
    start_time = time.time()
    minimised_gvalues = scaling_lbfgs(params, ndbins=15, nzbins=15, niter=2,
                                      n_absorption_positions=5, n_detector_bins=29,
                                      integration_method='summation_integration',
                                      decay_correction_rescaling=True,
                                      correction_order=['absorption', 'modulation', 'decay'])
    
                                
    total_time = time.time() - start_time
    print "algorithm execution time is " + str(total_time)

    '''save data and plot g-factors'''
    filename = filepath.strip('.pickle')+str('_scaled.txt')
    save_data(minimised_gvalues, filename)

    Rmeas = R_meas(minimised_gvalues.data_manager)
    Rpim = R_pim(minimised_gvalues.data_manager)
    print "R_meas of the (unmerged) data is %s" % (Rmeas)
    print "R_pim of the merged data is %s" % (Rpim)

    plot_data(minimised_gvalues.data_manager)
    plot_data_absorption(minimised_gvalues.data_manager)
    plot_data_modulation(minimised_gvalues.data_manager)
    
