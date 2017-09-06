from dials.array_family import flex
from dials.util.options import OptionParser
import minimiser_functions as mf
import data_manager_functions as dmf
from data_plotter import (plot_data, save_data, plot_data_absorption, 
                          plot_data_modulation)
from data_quality_assessment import R_meas, R_pim
import matplotlib.pyplot as plt
import time

def scaling_lbfgs(inputparams, gridding_parameters, niter=2,
                  integration_method='summation_integration',
                  decay_correction_rescaling=False,
                  correction_order=['absorption','modulation','decay']):
    """This algorithm loads a reflection table and json file and
    performs a Kabsch-scaling parameterisation using an lbfgs minimiser
    """
    if integration_method == 'profile_fitting':
        int_method = (str('intensity.prf.value'), str('intensity.prf.variance'))
    elif integration_method == 'summation_integration':
        int_method = (str('intensity.sum.value'), str('intensity.sum.variance'))
    else:
        raise ValueError('Incorrect choice of integration_method')

    loaded_reflections = dmf.Kabsch_Data_Manager(inputparams, int_method, gridding_parameters)

    """Filter out zero/negative values of d, variance, etc"""
    loaded_reflections.filter_data(int_method[1], -10.0, 0.0)
    loaded_reflections.filter_data('d', -1.0, 0.0)

    '''Map the indices to the asu and also sort the reflection table by miller index'''
    loaded_reflections.map_indices_to_asu()

    '''determine gridding index for scale parameters '''
    loaded_reflections.initialise_scale_factors()

    '''assign a unique reflection index to each group of reflections'''
    loaded_reflections.assign_h_index()
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
        for correction in correction_order:
            minimised = mf.optimiser(minimised.data_manager, correction=correction)

    #minimised.sort_data_for_output()
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
    gridding_parameters = {'ndbins':15,'nzbins':15,'n_absorption_positions':5,
                           'n_detector_bins':29}
    minimised_gvalues = scaling_lbfgs(params, gridding_parameters,  niter=2,
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

    #plot_data(minimised_gvalues.data_manager)
    #plot_data_absorption(minimised_gvalues.data_manager)
    #plot_data_modulation(minimised_gvalues.data_manager)
    
