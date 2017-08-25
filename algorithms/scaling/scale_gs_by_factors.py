from data_plotter import load_data
from data_quality_assessment import R_meas, R_pim
from dials.array_family import flex
from math import *
import matplotlib.pyplot as plt
import numpy as np

filename = "/Users/whi10850/Documents/test_data/integrate/13_integrated_scaled.txt"
data = load_data(filename)

data.sigma = 10000.0

print data.residual()
print list(data.x)[0:10]

def plot_data(data):
    y_label = 'd bin boundaries'
    y_ticks = ['%.3f' % x for x in data.data_manager.bin_boundaries['d']]
    x_label = 'z bin boundaries'
    x_ticks = data.data_manager.bin_boundaries['z_value'][::2]
    x_ticks = ['%.0f' % x for x in x_ticks]

    ndbins = len(data.data_manager.bin_boundaries['d'])-1
    nzbins = len(data.data_manager.bin_boundaries['z_value'])-1

    '''generate a plot of the result'''
    G_fin =  list(data.x)
    G_fin_2d=np.reshape(G_fin,(nzbins,ndbins)).T

    plt.figure(1)
    im=plt.imshow(G_fin_2d,cmap='viridis',origin='lower')
    plt.colorbar(im)
    plt.ylabel('d-value')
    plt.xlabel('time (z)')
    plt.yticks(np.arange(-0.5,ndbins),y_ticks)
    plt.xticks(np.arange(-0.5,nzbins,2),x_ticks)
    plt.title('$G_l$ correction factors using Kabsch method')
    plt.show()


def scale_by_Brel(data, Brel=0.7):
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
    print data.residual()
    plot_data(data)

def scale_by_constant_mult(data, m=0.15):

    data.data_manager.g_values = data.data_manager.g_values * m
    data.data_manager.calc_Ih()
    data.x = data.data_manager.g_values
    print list(data.x)[0:10]
    print data.residual()
    plot_data(data)

def scale_by_constant(data, c=2.0):
    nbins = data.data_manager.n_bins
    data.data_manager.g_values = data.data_manager.g_values + flex.double([c]*nbins)
    data.data_manager.calc_Ih()
    data.x = data.data_manager.g_values
    print list(data.x)[0:10]
    print data.residual()
    plot_data(data)



scale_by_Brel(data)
scale_by_constant_mult(data)
#scale_by_constant(data)

exit()


Rmeas = R_meas(data.data_manager)
Rpim = R_pim(data.data_manager)
print "R_meas of the (unmerged) data is %s" % (Rmeas)
print "R_pim of the merged data is %s" % (Rpim)

