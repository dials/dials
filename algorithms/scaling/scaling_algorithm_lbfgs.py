from dials.array_family import flex
from dials.util.options import (Importer, OptionParser)
import minimiser_functions as mf
import data_manager_functions as dmf
from data_plotter import plot_data, save_data


'''load data'''
filepath = "/Users/whi10850/Documents/test_data/integrate/13_integrated.pickle"
json_filepath = "/Users/whi10850/Documents/test_data/integrate/13_integrated_experiments.json"
parsestring = OptionParser(read_experiments = True, read_reflections = True,
                           check_format = False)
params, options = parsestring.parse_args([json_filepath,filepath])


'''Load reflection data in Data_Manager object'''
loaded_reflections = dmf.Data_Manager(params)

'''Filter out zero values of d, variance, etc'''
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

ndbins = 20; nzbins = 20
nlbins = ndbins * nzbins
resmax = (1.0 / (dmin**2))
resmin = (1.0 / (dmax**2))
resolution_bins = (flex.double(range(0, ndbins + 1)) * ((resmax - resmin)/ndbins) 
                + flex.double([resmin] * (ndbins + 1)))
d_bins = []
for n in resolution_bins:
    d_bins.append(1/(n**0.5))
d_bins = flex.double(d_bins[::-1])
z_bins = flex.double(range(0, nzbins + 1)) * zmax/nzbins

'''determine gridding index for scale parameters '''
loaded_reflections.bin_reflections(bin_parameters=[["d", d_bins],
                                                   ["z_value", z_bins]])


'''assign a unique reflection index to each group of reflections'''
loaded_reflections.assign_h_index()

'''Give initial scale values of 1.0'''
G0_values = flex.double([1.0] * (ndbins * nzbins))
loaded_reflections.set_g_values(G0_values)

'''call the optimiser on the Data Manager object'''
minimised = mf.optimiser(loaded_reflections, sigma = 0.05)		


filename = "/Users/whi10850/Documents/test_data/integrate/13_integrated_scaled.txt"
save_data(minimised, filename)

plot_data(datafile=filename)
