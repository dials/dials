from dials.array_family import flex
from dials.util.options import OptionParser
import minimiser_functions as mf
import data_manager_functions as dmf
from data_plotter import plot_data, save_data, plot_data_absorption, plot_data_modulation
from data_quality_assessment import R_meas, R_pim
import matplotlib.pyplot as plt
import time

def scaling_lbfgs(inputparams, ndbins, nzbins, sigma, niter=2):
    """This algorithm loads a reflection table and json file and
    performs a Kabsch-scaling parameterisation using an lbfgs minimiser
    """
    loaded_reflections = dmf.Data_Manager(inputparams)

    """Filter out zero values of d, variance, etc"""
    loaded_reflections.filter_data('intensity.sum.variance', -1.0, 0.0)
    loaded_reflections.filter_data('d', -1.0, 0.0)

    '''Map the indices to the asu and also sort the reflection table by miller index'''
    loaded_reflections.map_indices_to_asu()

    '''determine gridding index for scale parameters '''
    loaded_reflections.bin_reflections_dz(ndbins=ndbins, nzbins=nzbins)
    loaded_reflections.bin_reflections_absorption(npos=5)
    loaded_reflections.bin_reflections_modulation(ngridpoints=29)
    

    '''assign a unique reflection index to each group of reflections'''
    loaded_reflections.assign_h_index()
    loaded_reflections.create_index_tuple()
    loaded_reflections.scale_by_LP_and_dqe()

    '''call the optimiser on the Data Manager object''' 
    minimised = mf.optimiser(loaded_reflections, sigma, correction='absorption')
    minimised = mf.optimiser(minimised.data_manager, sigma, correction='modulation')
    minimised = mf.optimiser(minimised.data_manager, sigma, correction='decay')
    '''do invariant rescaling to physical B value'''
    minimised.data_manager.scale_gvalues()

    '''repeat n times - replace with convergence criteria and max iter instead?'''
    for i in range(0,niter):
        minimised = mf.optimiser(loaded_reflections, sigma, correction='absorption')
        minimised = mf.optimiser(minimised.data_manager, sigma, correction='modulation')
        minimised = mf.optimiser(minimised.data_manager, sigma, correction='decay')
        minimised.data_manager.scale_gvalues()
        
    return minimised


if __name__ == "__main__":
    filepath = "/Users/whi10850/Documents/test_data/integrate/13_integrated.pickle"
    json_filepath = "/Users/whi10850/Documents/test_data/integrate/13_integrated_experiments.json"
    parsestring = OptionParser(read_experiments=True, read_reflections=True,
                               check_format=False)
    params, options = parsestring.parse_args([json_filepath, filepath])

    '''do minimisation of g-factors'''
    start_time = time.time()
    minimised_gvalues = scaling_lbfgs(params, ndbins=15, nzbins=15, sigma=10000.0, niter=2)
    total_time = time.time() - start_time
    print "algorithm time is " + str(total_time)

    '''save data and plot g-factors'''
    filename = "/Users/whi10850/Documents/test_data/integrate/13_integrated_scaled.txt"
    save_data(minimised_gvalues, filename)

    Rmeas = R_meas(minimised_gvalues.data_manager)
    Rpim = R_pim(minimised_gvalues.data_manager)
    print "R_meas of the (unmerged) data is %s" % (Rmeas)
    print "R_pim of the merged data is %s" % (Rpim)

    plot_data(minimised_gvalues.data_manager)
    plot_data_absorption(minimised_gvalues.data_manager)
    plot_data_modulation(minimised_gvalues.data_manager)
    
