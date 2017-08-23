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
    G_fin =  list(data.x)

    y_label = 'd bin boundaries'
    y_ticks = ['%.3f' % x for x in data.data_manager.bin_boundaries['d']]
    x_label = 'z bin boundaries'
    x_ticks = data.data_manager.bin_boundaries['z_value'][::2]
    x_ticks = ['%.0f' % x for x in x_ticks]

    ndbins = len(data.data_manager.bin_boundaries['d'])-1
    nzbins = len(data.data_manager.bin_boundaries['z_value'])-1

    def plot_d_slices(G_fin):
        for i in range(0,ndbins):
            G_cut = G_fin[i:(nzbins*(ndbins-1))+i:ndbins]
            plt.plot(G_cut)
        plt.show()

    def plot_z_slices(G_fin):
        for i in range(0,nzbins):
            G_cut = G_fin[(ndbins*i):((i+1)*(ndbins))]
            plt.plot(G_cut)
        plt.show()

    '''generate a plot of the result'''

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

    #print 'Intensity, G_l, h, I_h'
    #data = load_data(filename)
    #for n, val in enumerate(data.data_manager.sorted_reflections['l_bin_index']):
    #    h = data.data_manager.sorted_reflections['h_index'][n]
    #    if h == 4554:
    #        h = data.data_manager.sorted_reflections['h_index'][n]
    #        print data.data_manager.sorted_reflections['intensity.sum.value'][n], data.x[val], val, data.data_manager.Ih_array[h]
    #exit()

    plot_data(datafile=filename)
