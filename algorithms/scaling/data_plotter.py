import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np

def save_data(minimised,filename):
    data_file = open(filename,'w')
    pickle.dump(minimised, data_file)
    data_file.close()

def load_data(filename):
    data_file = open(filename)
    data = pickle.load(data_file)
    data_file.close()
    return data

def plot_data(data_man):
    "takes in a data manager object"
    y_ticks = ['%.3f' % x for x in data_man.bin_boundaries['d']]
    x_ticks = data_man.bin_boundaries['z_value'][::2]
    x_ticks = ['%.0f' % x for x in x_ticks]

    ndbins = len(data_man.bin_boundaries['d'])-1
    nzbins = len(data_man.bin_boundaries['z_value'])-1
    '''generate a plot of the result'''
    G_fin = list(data_man.g_values)
    G_fin_2d = np.reshape(G_fin, (nzbins, ndbins)).T

    plt.figure(1)
    im = plt.imshow(G_fin_2d, cmap='viridis', origin='lower')
    plt.colorbar(im)
    plt.ylabel('d-value')
    plt.xlabel('time (z)')
    plt.yticks(np.arange(-0.5, ndbins), y_ticks)
    plt.xticks(np.arange(-0.5, nzbins, 2), x_ticks)
    plt.title('$G_l$ correction factors using Kabsch method')
    plt.savefig('Scaling_output_figure_lbfgs.png')
    plt.show()

def plot_data_absorption(data_man):
    "takes in a data manager object"
    #y_ticks = ['%.3f' % x for x in data_man.bin_boundaries['d']]
    x_ticks = data_man.bin_boundaries['z_value'][::2]
    x_ticks = ['%.0f' % x for x in x_ticks]

    
    nzbins = len(data_man.bin_boundaries['z_value'])-1
    nabsbins = len(data_man.g2_values)//nzbins
    '''generate a plot of the result'''
    G_fin = list(data_man.g2_values)
    G_fin_2d = np.reshape(G_fin, (nzbins, nabsbins)).T

    plt.figure(1)
    im = plt.imshow(G_fin_2d, cmap='viridis', origin='lower')
    plt.colorbar(im)
    plt.ylabel('detector position')
    plt.xlabel('time (z)')
    #plt.yticks(np.arange(-0.5, ndbins), y_ticks)
    plt.xticks(np.arange(-0.5, nzbins, 2), x_ticks)
    plt.title('$G_l$ correction factors using Kabsch method - \n detector position vs time')
    #plt.savefig('Scaling_output_figure_lbfgs.png')
    plt.show()

if __name__ == "__main__":
    datafile="/Users/whi10850/Documents/test_data/integrate/13_integrated_scaled.txt"
    data = load_data(filename = datafile)
    plot_data(data)
