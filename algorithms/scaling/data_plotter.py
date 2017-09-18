import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np
from data_quality_assessment import R_meas, R_pim

def save_data(minimised,filename):
  data_file = open(filename,'w')
  pickle.dump(minimised, data_file)
  data_file.close()

def load_data(filename):
  data_file = open(filename)
  data = pickle.load(data_file)
  data_file.close()
  return data

def plot_data_decay(data_man):
  "takes in a data manager object"
  y_ticks = ['%.3f' % x for x in data_man.bin_boundaries['d']]
  x_ticks = data_man.bin_boundaries['z_value'][::2]
  x_ticks = ['%.0f' % x for x in x_ticks]

  ndbins = len(data_man.bin_boundaries['d'])-1
  nzbins = len(data_man.bin_boundaries['z_value'])-1
  '''generate a plot of the result'''
  G_fin = list(data_man.g_decay)
  G_fin_2d = np.reshape(G_fin, (nzbins, ndbins)).T

  plt.figure(figsize=(7,7))
  im = plt.imshow(G_fin_2d, cmap='viridis', origin='lower')
  plt.colorbar(im)
  plt.ylabel('d-value')
  plt.xlabel('time (z)')
  plt.yticks(np.arange(-0.5, ndbins), y_ticks)
  plt.xticks(np.arange(-0.5, nzbins, 2), x_ticks)
  plt.title('Inverse scale factors for absorption correction')
  plt.savefig('g_resolution.png')
  #plt.show()

def plot_data_absorption(data_man):
  "takes in a data manager object"
  x_ticks = data_man.bin_boundaries['z_value'][::2]
  x_ticks = ['%.0f' % x for x in x_ticks]
  nzbins = len(data_man.bin_boundaries['z_value'])-1
  nabsbins = len(data_man.g_absorption)//nzbins
  '''generate a plot of the result'''
  G_fin = list(data_man.g_absorption)
  G_fin_2d = np.reshape(G_fin, (nzbins, nabsbins)).T

  plt.figure(figsize=(7,7))
  im = plt.imshow(G_fin_2d, cmap='viridis', origin='lower')
  plt.colorbar(im)
  plt.ylabel('detector position')
  plt.xlabel('time (z)')
  #plt.yticks(np.arange(-0.5, nabsbins), y_ticks)
  #plt.xticks(np.arange(-0.5, nzbins, 2), x_ticks)
  plt.title('Inverse scale factors for absorption correction')
  plt.savefig('g_absorption.png')
  #plt.show()

def plot_data_modulation(data_man):
  "takes in a data manager object"
  x_ticks = data_man.bin_boundaries['z_value'][::2]
  x_ticks = ['%.0f' % x for x in x_ticks]
  nbins = data_man.binning_parameters['n_detector_bins']
  
  '''generate a plot of the result'''
  G_fin = list(data_man.g_modulation)
  G_fin_2d = np.reshape(G_fin, (nbins, nbins))

  plt.figure(figsize=(7,7))
  im = plt.imshow(G_fin_2d, cmap='viridis', origin='lower')
  plt.colorbar(im)
  plt.ylabel('y')
  plt.xlabel('x')
  #plt.yticks(np.arange(-0.5, nbins), x_ticks)
  #plt.xticks(np.arange(-0.5, nbins), x_ticks)
  plt.title('Inverse scale factors for modulation correction')
  plt.savefig('g_modulation.png')
  #plt.show()

def plot_correction_at_detector_area(data_man, position):
  G_fin = list(data_man.g_absorption)
  npos = data_man.binning_parameters['n_absorption_positions']
  G_slice = G_fin[position::npos**2]
  plt.figure(1)
  plt.plot(G_slice)
  plt.xlabel('time (z)')
  plt.ylabel('$G_l$')
  plt.title('$G_l$ absorption correction factors as a function of time \n at the 13th detector position')
  plt.savefig('Scaling_output_figure_lbfgs_absorptionvstime.png')
  plt.show()

def plot_correction_at_resolution(data_man, position):
  G_fin = list(data_man.g_decay)
  nzbins = data_man.binning_parameters['n_z_bins']
  G_slice = G_fin[position::nzbins]
  plt.figure(1)
  plt.plot(G_slice)
  plt.xlabel('time (z)')
  plt.show()

def plot_absorption_correction_at_zbin(data_man, position):
  G_fin = list(data_man.g_absorption)
  npos = data_man.binning_parameters['n_absorption_positions']
  nzbins = data_man.binning_parameters['n_z_bins']
  G_slice = G_fin[position*npos*npos:(position+1)*npos*npos]
  G_slice_2d = np.reshape(G_slice, (npos, npos)).T
  plt.figure(1)
  im = plt.imshow(G_slice_2d, cmap='viridis', origin='lower')
  plt.colorbar(im)
  plt.xlabel('x bin')
  plt.ylabel('y bin')
  plt.title('$G_l$ correction factors for a $z$ bin')
  plt.savefig('Scaling_output_figure_lbfgs_detector.png')
  plt.show()

if __name__ == "__main__":
  datafile="/Users/whi10850/Documents/dials_scratch/jbe/scaling_code/test_data/x4_wide_integrated_scaled.pickle"
  data_man = load_data(filename = datafile)
  #data_man = data.data_manager
  
  #Rmeas = R_meas(data_man)
  #print "R_meas of the (unmerged) data is %s" % (Rmeas)
  #Rpim = R_pim(data_man)
  #print "R_pim of the merged data is %s" % (Rpim)
  plot_data(data_man)
  plot_data_absorption(data_man)
  plot_data_modulation(data_man)
  plot_absorption_correction_at_zbin(data_man, position=6)
  plot_correction_at_resolution(data_man, position=5)
  plot_correction_at_detector_area(data_man, position=13)
