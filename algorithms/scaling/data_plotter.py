import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np
from data_quality_assessment import R_meas, R_pim
import data_manager_functions as dmf
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

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
  y_ticks = data_man.bin_boundaries['d'][::2]
  y_ticks = ['%.2f' % x for x in y_ticks]
  x_ticks = data_man.bin_boundaries['z_value'][::2]
  x_ticks = ['%.0f' % x for x in x_ticks]
  ndbins = len(data_man.bin_boundaries['d'])-1
  nzbins = len(data_man.bin_boundaries['z_value'])-1
  '''generate a plot of the result'''
  G_fin = list(data_man.g_decay)
  G_fin_2d = np.reshape(G_fin, (nzbins, ndbins)).T


  plt.figure(figsize=(11,6))
  gs = gridspec.GridSpec(1, 2)
  ax1 = plt.subplot(gs[0, 0])
  ax2 = plt.subplot(gs[0, 1])
  ax2.hist(G_fin, 40, log=False)
  ax2.set_xlabel('Inverse scale factor')
  ax2.set_ylabel('Occurences')
  
  im = ax1.imshow(G_fin_2d, cmap='viridis', origin='lower')
  divider = make_axes_locatable(ax1)
  cax1 = divider.append_axes("right", size="5%", pad=0.05)
  cbar = plt.colorbar(im, cax=cax1)
  ax1.set_ylabel('d-value')
  ax1.set_xlabel('time (z)')
  ax1.set_yticks(np.arange(-0.5, ndbins, 2))
  ax1.set_yticklabels(y_ticks)
  ax1.set_xticks(np.arange(-0.5, nzbins, 2))
  ax1.set_xticklabels(x_ticks)
  ax1.set_title('Inverse scale factors for decay correction', fontsize=12)
  plt.tight_layout()
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

  plt.figure(figsize=(11,6))
  gs = gridspec.GridSpec(1, 2)
  ax1 = plt.subplot(gs[0, 0])
  ax2 = plt.subplot(gs[0, 1])
  ax2.hist(G_fin, 40, log=False)
  ax2.set_xlabel('Inverse scale factor')
  ax2.set_ylabel('Occurences')
  
  im = ax1.imshow(G_fin_2d, cmap='viridis', origin='lower')
  divider = make_axes_locatable(ax1)
  cax1 = divider.append_axes("right", size="5%", pad=0.05)
  cbar = plt.colorbar(im, cax=cax1)
  ax1.set_ylabel('detector position')
  ax1.set_xlabel('time (z)')
  #ax1.yticks(np.arange(-0.5, nabsbins), y_ticks)
  ax1.set_xticks(np.arange(-0.5, nzbins, 2))
  ax1.set_xticklabels(x_ticks)
  ax1.set_title('Inverse scale factors for absorption correction', fontsize=12)
  plt.tight_layout()
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

  plt.figure(figsize=(11,6))
  gs = gridspec.GridSpec(1, 2)
  ax1 = plt.subplot(gs[0, 0])
  ax2 = plt.subplot(gs[0, 1])
  ax2.hist(G_fin, 40, log=False)
  ax2.set_xlabel('Inverse scale factor')
  ax2.set_ylabel('Occurences')
  
  im = ax1.imshow(G_fin_2d, cmap='viridis', origin='lower')
  divider = make_axes_locatable(ax1)
  cax1 = divider.append_axes("right", size="5%", pad=0.05)
  cbar = plt.colorbar(im, cax=cax1)
  ax1.set_ylabel('y')
  ax1.set_xlabel('x')
  #plt.yticks(np.arange(-0.5, nbins), x_ticks)
  #plt.xticks(np.arange(-0.5, nbins), x_ticks)
  ax1.set_title('Inverse scale factors for modulation correction', fontsize=12)
  plt.tight_layout()
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
  
  plot_data_decay(data_man)
  plot_data_absorption(data_man)
  plot_data_modulation(data_man)
  #plot_absorption_correction_at_zbin(data_man, position=6)
  #plot_correction_at_resolution(data_man, position=5)
  #plot_correction_at_detector_area(data_man, position=13)
