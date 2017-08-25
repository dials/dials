from data_plotter import load_data
from data_quality_assessment import R_meas, R_pim
from dials.array_family import flex
from math import *

filename = "/Users/whi10850/Documents/test_data/integrate/13_integrated_scaled.txt"
data = load_data(filename)

data.sigma = 10000.0

print data.residual()
print list(data.x)[0:10]


def scale_by_Brel(data, Brel=-1.0):
    d_bin_boundaries = data.data_manager.bin_boundaries['d']
    d_values = flex.double([])
    for i in range(0,len(d_bin_boundaries)-1):
        d_values.append((d_bin_boundaries[i] + d_bin_boundaries[i+1])/2.0)

    scaling_factors = []
    for i in range(0,data.data_manager.nzbins):
        scaling_factors += flex.exp(Brel/(d_values**2))

    scaling_factors = flex.double(scaling_factors)

    data.data_manager.g_values = data.data_manager.g_values * scaling_factors
    data.data_manager.calc_Ih()
    data.x = data.data_manager.g_values
    print list(data.x)[0:10]
    print data.residual()

def scale_by_constant_mult(data, m=2.0):

    data.data_manager.g_values = data.data_manager.g_values * m
    data.data_manager.calc_Ih()
    data.x = data.data_manager.g_values
    print list(data.x)[0:10]
    print data.residual()

def scale_by_constant(data, c=2.0):
    nbins = data.data_manager.n_bins
    data.data_manager.g_values = data.data_manager.g_values + flex.double([c]*nbins)
    data.data_manager.calc_Ih()
    data.x = data.data_manager.g_values
    print list(data.x)[0:10]
    print data.residual()



scale_by_Brel(data)
scale_by_constant_mult(data)
scale_by_constant(data)

exit()


Rmeas = R_meas(data.data_manager)
Rpim = R_pim(data.data_manager)
print "R_meas of the (unmerged) data is %s" % (Rmeas)
print "R_pim of the merged data is %s" % (Rpim)

