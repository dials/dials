from dials.array_family import flex
from dials.util.options import OptionParser
import minimiser_functions as mf
import data_manager_functions as dmf
from data_plotter import plot_data, save_data
from data_quality_assessment import R_meas, R_pim

def scaling_lbfgs(inputparams, ndbins, nzbins, sigma):
    """This algorithm loads a reflection table and json file and
    performs a Kabsch-scaling parameterisation using an lbfgs minimiser
    """

    loaded_reflections = dmf.Data_Manager(inputparams)

    """Filter out zero values of d, variance, etc"""
    loaded_reflections.filter_data('intensity.sum.variance', -1.0, 0.0)
    loaded_reflections.filter_data('d', -1.0, 0.0)

    '''Map the indices to the asu and also sort the reflection table by miller index'''
    loaded_reflections.map_indices_to_asu()

    '''Determine binning/gridding boundaries'''
    zmax = max(loaded_reflections.filtered_reflections['z_value'])
    zmin = min(loaded_reflections.filtered_reflections['z_value'])
    dmax = max(loaded_reflections.filtered_reflections['d'])
    dmin = min(loaded_reflections.filtered_reflections['d'])
    print 'dmin,dmax,zmin,zmax = %f, %f, %f, %f' % (dmin, dmax, zmin, zmax)

    nlbins = ndbins * nzbins
    resmax = (1.0 / (dmin**2))
    resmin = (1.0 / (dmax**2))
    resolution_bins = (flex.double(range(0, ndbins + 1)) * ((resmax - resmin)/ndbins)
                       + flex.double([resmin] * (ndbins + 1)))
    d_bins = []
    for bin_boundary in resolution_bins:
        d_bins.append(1/(bin_boundary**0.5))
    d_bins = flex.double(d_bins[::-1])
    z_bins = flex.double(range(0, nzbins + 1)) * zmax/nzbins

    '''determine gridding index for scale parameters '''
    loaded_reflections.bin_reflections(bin_parameters=[["d", d_bins],
                                                       ["z_value", z_bins]])

    '''assign a unique reflection index to each group of reflections'''
    loaded_reflections.assign_h_index()

    '''Give initial scale values of 1.0'''
    initial_g_values = flex.double([1.0] * (nlbins))
    loaded_reflections.set_g_values(initial_g_values)

    '''call the optimiser on the Data Manager object'''
    minimised = mf.optimiser(loaded_reflections, sigma)

    return minimised

if __name__ == "__main__":
    filepath = "/Users/whi10850/Documents/test_data/integrate/13_integrated.pickle"
    json_filepath = "/Users/whi10850/Documents/test_data/integrate/13_integrated_experiments.json"
    parsestring = OptionParser(read_experiments=True, read_reflections=True,
                               check_format=False)
    params, options = parsestring.parse_args([json_filepath, filepath])

    '''do minimisation of g-factors'''
    minimised_gvalues = scaling_lbfgs(params, ndbins=20, nzbins=20, sigma=0.05)

    '''save data and plot g-factors'''
    filename = "/Users/whi10850/Documents/test_data/integrate/13_integrated_scaled.txt"
    save_data(minimised_gvalues, filename)

    Rmeas = R_meas(minimised_gvalues.data_manager)
    Rpim = R_pim(minimised_gvalues.data_manager)
    print "R_meas of the (unmerged) data is %s" % (Rmeas)
    print "R_pim of the merged data is %s" % (Rpim)


    plot_data(datafile=filename)
