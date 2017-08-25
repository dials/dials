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

def plot_data(datafile):
    data = load_data(filename = datafile)
    

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
    plt.savefig('Scaling_output_figure_lbfgs.png')
    plt.figure(2)
    plt.plot(data.residuals)
    plt.ylabel('Residual')
    plt.show()

if __name__ == "__main__":
    filename="/Users/whi10850/Documents/test_data/integrate/13_integrated_scaled.txt"

    plot_data(datafile=filename)
