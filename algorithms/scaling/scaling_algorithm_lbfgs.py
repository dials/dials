import time
import sys
from dials.array_family import flex
from dials.util.options import OptionParser
import minimiser_functions as mf
import data_manager_functions as dmf
from data_plotter import (plot_data, save_data, plot_data_absorption,
                          plot_data_modulation)
from data_quality_assessment import R_meas, R_pim
import matplotlib.pyplot as plt


def scaling_lbfgs(inputparams, gridding_parameters, scaling_options):
    """This algorithm loads a reflection table and json file and
    performs a Kabsch-scaling parameterisation using an lbfgs minimiser
    """
    '''handling of choice of integration method'''
    if scaling_options['integration_method'] == 'profile_fitting':
        int_method = (str('intensity.prf.value'), str('intensity.prf.variance'))
    elif scaling_options['integration_method'] == 'summation_integration':
        int_method = (str('intensity.sum.value'), str('intensity.sum.variance'))
    else:
        raise ValueError('Incorrect choice of integration_method')

    '''create a data manager object'''
    loaded_reflections = dmf.Kabsch_Data_Manager(inputparams, int_method,
                                                 gridding_parameters)

    """Filter out zero/negative values of d, variance, etc"""
    loaded_reflections.filter_data(int_method[1], -10.0, 0.0)
    loaded_reflections.filter_data('d', -1.0, 0.0)
    #loaded_reflections.filter_I_sigma(1.0)

    '''Map the indices to the asu and also sort the reflection table by miller index'''
    loaded_reflections.map_indices_to_asu()

    '''determine gridding index for scale parameters '''
    loaded_reflections.initialise_scale_factors()

    '''assign a unique reflection index to each group of reflections'''
    loaded_reflections.assign_h_index()
    loaded_reflections.scale_by_LP_and_dqe()

    '''call the optimiser on the Data Manager object'''
    correction_order = scaling_options['correction_order']
    decay_correction_rescaling = scaling_options['decay_correction_rescaling']
    minimised = mf.optimiser(loaded_reflections, correction=correction_order[0],
                             decay_correction_rescaling=decay_correction_rescaling)
    minimised = mf.optimiser(minimised.data_manager, correction=correction_order[1],
                             decay_correction_rescaling=decay_correction_rescaling)
    minimised = mf.optimiser(minimised.data_manager, correction=correction_order[2],
                             decay_correction_rescaling=decay_correction_rescaling)
    
    #minimised.data_manager.reject_outliers(tolerance=6.0,niter=4)
    #minimised.data_manager.assign_h_index()
    '''repeat n times - replace with convergence criteria and max iter instead?'''
    for _ in range(0, scaling_options['n_cycles'] - 1):
        for correction in correction_order:
            minimised = mf.optimiser(minimised.data_manager, correction=correction)
        #minimised.data_manager.reject_outliers(tolerance=6.0,niter=1)
        #minimised.data_manager.assign_h_index()
    minimised.data_manager.sorted_reflections['inverse_scale_factor'] = minimised.data_manager.scale_factors
    
    return minimised



if __name__ == "__main__":
    '''load datafiles - can be replaced with sys.argv commands'''
    parsestring = OptionParser(read_experiments=True, read_reflections=True,
                               check_format=False)
    #filepath = "/Users/whi10850/Documents/test_data/integrate/13_integrated.pickle"
    #json_filepath = "/Users/whi10850/Documents/test_data/integrate/13_integrated_experiments.json"
    #params, options = parsestring.parse_args([json_filepath, filepath])
    '''hack to 'overload' the phil_scope methods and extract my own parameters
    into a dict'''
    scaling_options_dict = {}
    if len(sys.argv) > 3:
        try:
            input_parameters = sys.argv[3:]
            for i in input_parameters:
                scaling_options_dict[i.split('=')[0]] = int(i.split('=')[1])
            print "input parameters read in: " +str(scaling_options_dict)
        except:
            pass
    #expects an integrated.pickle and an experiments.json
    params, options = parsestring.parse_args(sys.argv[1:3])
    start_time = time.time()

    #default parameters
    gridding_parameters = {'ndbins':20, 'nzbins':20, 'n_detector_bins':37}
    scaling_options = {'n_cycles':1, 'integration_method':'profile_fitting',
                       'decay_correction_rescaling':True,
                       'correction_order':['absorption', 'decay', 'modulation']}
    #hacky overload of default parameters from command line input
    for parameter, value in scaling_options_dict.iteritems():
        if parameter in gridding_parameters:
            gridding_parameters[parameter] = value
        elif parameter in scaling_options:
            scaling_options[parameter] = value
        else:
            print "parameter " +str(parameter) +" not recognised, skipping assignment"

    '''do minimisation of g-factors'''
    minimised_optimiser = scaling_lbfgs(params, gridding_parameters, scaling_options)

    total_time = time.time() - start_time
    print "algorithm execution time is " + str(total_time)

    '''save data and plot g-factors'''
    #filename = sys.argv[1].strip('.pickle')+str('_scaled.pickle')
    #save_data(minimised_optimiser, filename)
    #minimised_optimiser.data_manager.clean_sorted_reflections()
    #save_data(minimised_optimiser.data_manager, filename)
    filename = sys.argv[1].strip('.pickle')+str('_scaled.pickle')
    #minimised_optimiser.data_manager.create_fake_dataset()
    minimised_optimiser.data_manager.save_sorted_reflections(filename)

    Rmeas = R_meas(minimised_optimiser.data_manager)
    Rpim = R_pim(minimised_optimiser.data_manager)
    print "R_meas is %s" % (Rmeas)
    print "R_pim is %s" % (Rpim)

    plot_data(minimised_optimiser.data_manager)
    plot_data_absorption(minimised_optimiser.data_manager)
    plot_data_modulation(minimised_optimiser.data_manager)
