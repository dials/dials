import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np
from data_quality_assessment import R_meas, R_pim
import data_manager_functions as dmf
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cctbx.array_family import flex

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
  ndbins = data_man.g_decay.n1_parameters
  nzbins = data_man.g_decay.n2_parameters
  '''create a grid of x and y points and use these to generate scale factors'''
  rel_values_1 = np.arange(0, int(max(data_man.sorted_reflections['normalised_res_values'])) + 1, 0.1)
  rel_values_2 = np.arange(0, int(max(data_man.sorted_reflections['normalised_time_values'])) + 1, 0.1)
  (n1, n2) = (len(rel_values_1), len(rel_values_2))
  rel_values_1 = np.tile(rel_values_1, n2)
  rel_values_2 = np.repeat(rel_values_2, n1)
  test_scale_factor = dmf.SmoothScaleFactor_2D(1.0, ndbins, nzbins)
  test_scale_factor.set_scale_factors(data_man.g_decay.get_scale_factors())
  test_scale_factor.set_normalised_values(rel_values_1, rel_values_2)
  scales = test_scale_factor.calculate_smooth_scales()
  G_fin_2d = np.reshape(list(scales), (n2, n1)).T
  G_fin_2d = G_fin_2d[1:-1,1:-1]

  plt.figure(figsize=(11,6))
  gs = gridspec.GridSpec(1, 2)
  ax1 = plt.subplot(gs[0, 0])
  ax2 = plt.subplot(gs[0, 1])
  ax2.hist(scales, 40, log=False)
  ax2.set_xlabel('Inverse scale factor')
  ax2.set_ylabel('Occurences')
  
  im = ax1.imshow(G_fin_2d, cmap='viridis', origin='lower')
  divider = make_axes_locatable(ax1)
  cax1 = divider.append_axes("right", size="5%", pad=0.05)
  cbar = plt.colorbar(im, cax=cax1)
  ax1.set_ylabel('d-value')
  ax1.set_xlabel('time (z)')
  #ax1.set_yticks(np.arange(-0.5, ndbins, 2))
  #ax1.set_yticklabels(y_ticks)
  #ax1.set_xticks(np.arange(-0.5, nzbins, 2))
  #ax1.set_xticklabels(x_ticks)
  ax1.set_title('Inverse scale factors for decay correction', fontsize=12)
  plt.tight_layout()
  plt.savefig('g_resolution.png')
  #plt.show()

def plot_data_absorption(data_man):
  "takes in a data manager object"
  x_ticks = data_man.bin_boundaries['z_value'][::2]
  x_ticks = ['%.0f' % x for x in x_ticks]
  nzbins = len(data_man.bin_boundaries['z_value'])-1
  nabsbins = len(data_man.g_absorption.get_scale_factors())//nzbins
  '''generate a plot of the result'''
  G_fin = list(data_man.g_absorption.get_scale_factors())
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

  nxbins = data_man.g_modulation.n1_parameters
  nybins = data_man.g_modulation.n2_parameters

  '''create a grid of x and y points and use these to generate scale factors'''
  rel_values_1 = np.arange(0, int(max(data_man.sorted_reflections['normalised_x_values'])) + 1, 0.2)
  rel_values_2 = np.arange(0, int(max(data_man.sorted_reflections['normalised_y_values'])) + 1, 0.2)
  (n1, n2) = (len(rel_values_1), len(rel_values_2))
  rel_values_1 = np.tile(rel_values_1, n2)
  rel_values_2 = np.repeat(rel_values_2, n1)
  test_scale_factor = dmf.SmoothScaleFactor_2D(1.0, nxbins, nybins)
  test_scale_factor.set_scale_factors(data_man.g_modulation.get_scale_factors())
  test_scale_factor.set_normalised_values(rel_values_1, rel_values_2)
  scales = test_scale_factor.calculate_smooth_scales()
  G_fin_2d = np.reshape(list(scales), (n2, n1))
  G_fin_2d = G_fin_2d[1:-1,1:-1]

  plt.figure(figsize=(11,6))
  gs = gridspec.GridSpec(1, 2)
  ax1 = plt.subplot(gs[0, 0])
  ax2 = plt.subplot(gs[0, 1])
  ax2.hist(scales, 40, log=False)
  ax2.set_xlabel('Inverse scale factor')
  ax2.set_ylabel('Occurences')
  
  im = ax1.imshow(G_fin_2d, cmap='viridis', origin='lower')
  divider = make_axes_locatable(ax1)
  cax1 = divider.append_axes("right", size="5%", pad=0.05)
  cbar = plt.colorbar(im, cax=cax1)
  ax1.set_ylabel('y')
  ax1.set_xlabel('x')
  ax1.set_title('Inverse scale factors for modulation correction', fontsize=12)
  plt.tight_layout()
  plt.savefig('g_modulation.png')

def plot_correction_at_multiple_detector_areas(data_man, positions):
  n = len(positions)
  plt.figure(figsize=(12, 8))
  gs = gridspec.GridSpec(2, 3)
  min_sf = min(data_man.g_absorption.get_scale_factors())
  max_sf = max(data_man.g_absorption.get_scale_factors())
  for i, position in enumerate(positions):
    ax = plt.subplot(gs[i//3, i%3])
    G = calc_correction_at_detector_area(data_man, position)
    im = ax.imshow(G, cmap='viridis', origin='lower',vmin=min_sf, vmax=max_sf)
    divider = make_axes_locatable(ax)
    cax1 = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax1)
    ax.set_ylabel('y')
    ax.set_xlabel('x')
    #ax.set_yticks(np.arange(data_man.g_absorption.ny_parameters*10, 10), 
    #  np.arange(data_man.g_absorption.ny_parameters))
    #ax.set_xticks(np.arange(data_man.g_absorption.nx_parameters*10, 10), 
    #  np.arange(data_man.g_absorption.nx_parameters))
    ax.set_title('Absorption correction factor surface at time bin %s' % (position), fontsize=7)
  plt.tight_layout()
  plt.savefig('g_absorption_surfaces.png')


def calc_correction_at_detector_area(data_man, position):
  nxbins = data_man.g_absorption.nx_parameters
  nybins = data_man.g_absorption.ny_parameters
  #form an xy grid
  rel_values_1 = np.arange(0, int(max(data_man.sorted_reflections['normalised_x_abs_values'])) + 1, 0.1)
  rel_values_2 = np.arange(0, int(max(data_man.sorted_reflections['normalised_y_abs_values'])) + 1, 0.1)
  (n1, n2) = (len(rel_values_1), len(rel_values_2))
  rel_values_1 = np.tile(rel_values_1, n2)
  rel_values_2 = np.repeat(rel_values_2, n1)
  rel_values_3 = flex.double([0.0]*len(rel_values_2))

  test_scale_factor = dmf.SmoothScaleFactor_GridAbsorption(1.0, nxbins, nybins, 1)
  test_scale_factor.set_scale_factors(data_man.g_absorption.get_scale_factors()[
    position*nxbins*nybins:(position+1)*nxbins*nybins])
  test_scale_factor.set_normalised_values(rel_values_1, rel_values_2, rel_values_3)
  scales = test_scale_factor.calculate_smooth_scales()
  G_fin_2d = np.reshape(list(scales), (n2, n1))

  return G_fin_2d


def plot_correction_at_resolution(data_man, position):
  G_fin = list(data_man.g_decay.get_scale_factors())
  nzbins = data_man.binning_parameters['n_z_bins']
  G_slice = G_fin[position::nzbins]
  plt.figure(1)
  plt.plot(G_slice)
  plt.xlabel('time (z)')
  plt.show()

def plot_absorption_correction_at_zbin(data_man, position):
  G_fin = list(data_man.g_absorption.get_scale_factors())
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

def plot_smooth_scales(data_man):                                           
  rel_values = np.arange(0, int(max(data_man.sorted_reflections['normalised_rotation_angle'])) + 1, 0.1)
  test_scale_factor = dmf.SmoothScaleFactor_1D(1.0, data_man.n_g_scale_params)
  test_scale_factor.set_scale_factors(data_man.g_scale.get_scale_factors())
  test_scale_factor.set_normalised_values(rel_values)
  scales = test_scale_factor.calculate_smooth_scales()
    
  plt.subplot(2,1,1)
  plt.title('Smooth scale factors')
  plt.plot(rel_values, scales)
  plt.ylabel('Scale term')
  plt.xlabel('Normalised rotation angle')

  rel_values = np.arange(0, int(max(data_man.sorted_reflections['normalised_time_values'])) + 1, 0.1)
  test_decay_factor = dmf.SmoothScaleFactor_1D(0.0, data_man.n_g_decay_params)
  test_decay_factor.set_scale_factors(data_man.g_decay.get_scale_factors())
  test_decay_factor.set_normalised_values(rel_values)
  B_rel_values = test_decay_factor.calculate_smooth_scales()
  plt.subplot(2,1,2)
  plt.ylabel('Relative B factor')
  plt.xlabel('Normalised time value')
  plt.plot(rel_values, B_rel_values)
  plt.tight_layout()
  plt.savefig('Smooth_scale_factors')

if __name__ == "__main__":
  datafile="/Users/whi10850/Documents/dials_scratch/jbe/scaling_code/test_data/x4_wide_integrated_scaled.pickle"
  data_man = load_data(filename = datafile)
  
  plot_data_decay(data_man)
  plot_data_absorption(data_man)
  plot_data_modulation(data_man)
