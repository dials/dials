import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np
import data_manager_functions as dmf
import scale_factor as SF
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
#from cctbx.array_family import flex
from dials.array_family import flex
import copy as copy

def save_data(minimised,filename):
  data_file = open(filename,'w')
  pickle.dump(minimised, data_file)
  data_file.close()

def load_data(filename):
  data_file = open(filename)
  data = pickle.load(data_file)
  data_file.close()
  return data

def plot_data_decay(data_man, outputfile=None):
  ndbins = data_man.g_decay.n1_parameters
  nzbins = data_man.g_decay.n2_parameters
  '''create a grid of x and y points and use these to generate scale factors'''
  max_res = int(max(data_man.reflection_table['normalised_res_values'])) + 1
  max_time = int(max(data_man.reflection_table['normalised_time_values'])) + 1
  rel_values_1 = np.arange(0, max_res+0.1, 0.1)
  rel_values_2 = np.arange(0, max_time+0.1, 0.1)
  (n1, n2) = (len(rel_values_1), len(rel_values_2))
  rel_values_1 = np.tile(rel_values_1, n2)
  rel_values_2 = np.repeat(rel_values_2, n1)
  test_scale_factor = SF.SmoothScaleFactor_2D(1.0, ndbins, nzbins)
  test_scale_factor.set_scale_factors(data_man.g_decay.get_scale_factors())
  test_scale_factor.set_normalised_values(rel_values_1, rel_values_2)
  test_scale_factor.calculate_scales_and_derivatives()
  scales = test_scale_factor.scales
  G_fin_2d = np.reshape(list(scales), (n2, n1)).T
  #G_fin_2d = G_fin_2d[1:-1,1:-1]
  data_man.g_decay.calculate_scales_and_derivatives()
  G_fin = list(data_man.g_decay.scales)

  plt.figure(figsize=(11,6))
  gs = gridspec.GridSpec(1, 2)
  ax1 = plt.subplot(gs[0, 0])
  ax2 = plt.subplot(gs[0, 1])
  ax2.hist(G_fin, 40, log=False)
  ax2.set_xlabel('Inverse scale factor')
  ax2.set_ylabel('Occurences in reflections_for_scaling')

  resmax = (1.0 / (min(data_man.reflection_table['d'])**2))
  resmin = (1.0 / (max(data_man.reflection_table['d'])**2))

  #print "dmin is %s" % (min(data_man.reflection_table['d']))
  #print "dmax is %s" % (max(data_man.reflection_table['d']))
  #determine boundaries in res
  resbin_boundaries = np.arange(resmin, resmax, 2*(resmax - resmin)/(ndbins-1))
  #print resbin_boundaries
  dbin_boundaries = 1.0/(resbin_boundaries**0.5)
  #print dbin_boundaries
  dbin_boundaries = ['%.3f' % x for x in dbin_boundaries]

  im = ax1.imshow(G_fin_2d, cmap='viridis', origin='upper', aspect='auto')
  divider = make_axes_locatable(ax1)
  cax1 = divider.append_axes("right", size="5%", pad=0.05)
  cbar = plt.colorbar(im, cax=cax1)
  ax1.set_ylabel('d-value')
  ax1.set_xlabel('Normalised time value')
  ax1.set_yticks(np.arange(0, (max_res * 10)+0.01, 20))
  ax1.set_yticklabels(dbin_boundaries)
  ax1.set_xticks(np.arange(0, (max_time * 10)+0.01, 20))
  ax1.set_xticklabels(np.arange(0, max_time+0.01, 2))
  ax1.set_title('Inverse scale factors for decay correction', fontsize=12)
  plt.tight_layout()
  if outputfile:
    plt.savefig(outputfile)
  else:
    plt.savefig('g_decay.png')
  #plt.show()

def plot_data_absorption(data_man, outputfile=None):
  "takes in a data manager object"

  nzbins = data_man.g_absorption.ntime_parameters
  nabsbins = data_man.g_absorption.nx_parameters * data_man.g_absorption.ny_parameters

  '''generate a plot of the result'''
  G_fin = list(data_man.g_absorption.get_scale_factors())
  G_fin_2d = np.reshape(G_fin, (nzbins, nabsbins)).T
  data_man.g_absorption.calculate_scales_and_derivatives()
  G_fin = list(data_man.g_absorption.scales)

  plt.figure(figsize=(11,6))
  gs = gridspec.GridSpec(1, 2)
  ax1 = plt.subplot(gs[0, 0])
  ax2 = plt.subplot(gs[0, 1])
  ax2.hist(G_fin, 40, log=False)
  ax2.set_xlabel('Inverse scale factor')
  ax2.set_ylabel('Occurences in reflections_for_scaling')

  im = ax1.imshow(G_fin_2d, cmap='viridis', origin='lower')
  divider = make_axes_locatable(ax1)
  cax1 = divider.append_axes("right", size="5%", pad=0.05)
  cbar = plt.colorbar(im, cax=cax1)
  ax1.set_ylabel('Detector position')
  ax1.set_xlabel('Normalised time value')
  
  #ax1.yticks(np.arange(-0.5, nabsbins), y_ticks)
  #ax1.set_xticks(np.arange(-0.5, nzbins, 2))
  #ax1.set_xticklabels(x_ticks)
  ax1.set_title('Inverse scale factors for absorption correction', fontsize=12)
  plt.tight_layout()
  if outputfile:
    plt.savefig(outputfile)
  else:
    plt.savefig('g_absorption.png')
  #plt.show()

def plot_data_modulation(data_man, outputfile=None):

  nxbins = data_man.g_modulation.n1_parameters
  nybins = data_man.g_modulation.n2_parameters
  max_x = int(max(data_man.reflection_table['normalised_x_values'])) + 1
  max_y = int(max(data_man.reflection_table['normalised_y_values'])) + 1
  '''create a grid of x and y points and use these to generate scale factors'''
  rel_values_1 = np.arange(0, max_x+0.2, 0.2)
  rel_values_2 = np.arange(0, max_y+0.2, 0.2)
  (n1, n2) = (len(rel_values_1), len(rel_values_2))
  rel_values_1 = np.tile(rel_values_1, n2)
  rel_values_2 = np.repeat(rel_values_2, n1)
  test_scale_factor = SF.SmoothScaleFactor_2D(1.0, nxbins, nybins)
  test_scale_factor.set_scale_factors(data_man.g_modulation.get_scale_factors())
  test_scale_factor.set_normalised_values(rel_values_1, rel_values_2)
  test_scale_factor.calculate_scales_and_derivatives()
  scales = test_scale_factor.scales
  G_fin_2d = np.reshape(list(scales), (n2, n1))

  plt.figure(figsize=(11,6))
  gs = gridspec.GridSpec(1, 2)
  ax1 = plt.subplot(gs[0, 0])
  ax2 = plt.subplot(gs[0, 1])
  ax2.hist(scales, 40, log=False)
  ax2.set_xlabel('Inverse scale factor')
  ax2.set_ylabel('Occurences in reflections_for_scaling')

  im = ax1.imshow(G_fin_2d, cmap='viridis', origin='lower')
  divider = make_axes_locatable(ax1)
  cax1 = divider.append_axes("right", size="5%", pad=0.05)
  cbar = plt.colorbar(im, cax=cax1)
  ax1.set_yticks(np.arange(0, (max_x * 5)+0.01, 10))
  ax1.set_yticklabels(np.arange(0, max_x+0.01, 2))
  ax1.set_xticks(np.arange(0, (max_y * 5)+0.01, 10))
  ax1.set_xticklabels(np.arange(0, max_y+0.01, 2))
  ax1.set_ylabel('Normalised y value')
  ax1.set_xlabel('Normalised x value')
  ax1.set_title('Inverse scale factors for modulation correction', fontsize=12)
  plt.tight_layout()
  if outputfile:
    plt.savefig(outputfile)
  else:
    plt.savefig('g_modulation.png')

def plot_correction_at_multiple_detector_areas(data_man, positions, outputfile=None):
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
    nx = data_man.g_absorption.nx_parameters
    ny = data_man.g_absorption.ny_parameters
    ax.set_yticks(np.arange(-0.5, (ny * 10) - 0.5, 10))
    ax.set_ylabel('Normalised y position')
    ax.set_yticklabels(np.arange(0, ny))
    ax.set_xticks(np.arange(-0.5, (nx * 10) - 0.5, 10))
    ax.set_xlabel('Normalised x position')
    ax.set_xticklabels(np.arange(0, nx))
    ax.set_title('Absorption correction factor surface at normalised time %s' % (float(position)+0.5), fontsize=7)
    #print "successfully plotted positon %s" % position
  plt.tight_layout()
  if outputfile:
    plt.savefig(outputfile)
  else:
    plt.savefig('g_absorption_surfaces.png')


def calc_correction_at_detector_area(data_man, position):
  nxbins = data_man.g_absorption.nx_parameters
  nybins = data_man.g_absorption.ny_parameters
  ntimebins = data_man.g_absorption.ntime_parameters
  #form an xyz grid, for two planes - going to plot at z = int + 0.5 between params
  rel_values_1 = np.arange(0, int(max(data_man.reflection_table['normalised_x_abs_values'])) + 1, 0.1)
  rel_values_2 = np.arange(0, int(max(data_man.reflection_table['normalised_y_abs_values'])) + 1, 0.1)
  (n1, n2) = (len(rel_values_1), len(rel_values_2))
  rel_values_1 = np.tile(rel_values_1, n2)
  rel_values_2 = np.repeat(rel_values_2, n1)
  rel_values_3 = flex.double([float(position)+0.5]*len(rel_values_2)) #why -0.5 - check.
  test_scale_factor = copy.deepcopy(data_man.g_absorption)
  test_scale_factor.set_normalised_values(rel_values_1, rel_values_2, rel_values_3)
  test_scale_factor.calculate_scales_and_derivatives()
  scales = test_scale_factor.scales
  G_fin_2d = np.reshape(list(scales), (n2, n1))
  return G_fin_2d


def plot_smooth_scales(data_man, outputfile=None):
  rel_values = np.arange(0, int(max(data_man.reflection_table['normalised_rotation_angle'])) + 1, 0.1)
  scale_SFs = data_man.g_scale.get_scale_factors()
  n_g_scale_params = len(scale_SFs)
  test_scale_factor = SF.SmoothScaleFactor_1D(1.0, n_g_scale_params)
  test_scale_factor.set_scale_factors(scale_SFs)
  test_scale_factor.set_normalised_values(rel_values)
  test_scale_factor.calculate_scales()
  scales = test_scale_factor.scales
  plt.figure(figsize=(14, 8))
  plt.subplot(2,1,1)
  plt.title('Smooth scale factors')
  plt.plot(rel_values, scales)
  plt.ylabel('Scale term')
  plt.xlabel('Normalised rotation angle')

  if data_man.params.parameterisation.decay_term:
    rel_values = np.arange(0, int(max(data_man.reflection_table['normalised_time_values'])) + 1, 0.1)
    decay_SFs = data_man.g_decay.get_scale_factors()
    n_g_decay_params = len(decay_SFs)
    test_decay_factor = SF.SmoothScaleFactor_1D(0.0, n_g_decay_params)
    test_decay_factor.Vr = 0.5 ##HACK - set to match that of SmoothScaleFactor_1D_Bfactor
    test_decay_factor.set_scale_factors(decay_SFs)
    test_decay_factor.set_normalised_values(rel_values)
    test_decay_factor.calculate_scales()
    B_rel_values = test_decay_factor.scales
    plt.subplot(2,1,2)
    plt.ylabel('Relative B factor (' + r'$\AA^{2}$'+')')
    plt.xlabel('Normalised time value')
    plt.plot(rel_values, B_rel_values)
  #plt.tight_layout()
  if outputfile:
    plt.savefig(outputfile)
  else:
    plt.savefig('smooth_scale_factors.png')

def plot_absorption_surface(data_man, outputfile=None):
  params = data_man.g_absorption.get_scale_factors()

  from scitbx import math
  from scitbx.array_family import flex
  import math as pymath
  order = data_man.params.parameterisation.lmax
  lfg =  math.log_factorial_generator(2 * order + 1)
  STEPS = 50
  phi = np.linspace(0, 2 * np.pi, 2*STEPS)
  theta = np.linspace(0, np.pi, STEPS)
  THETA, PHI = np.meshgrid(theta, phi)
  lmax = data_man.params.parameterisation.lmax
  Intensity = np.ones(THETA.shape)
  counter = 0
  sqrt2 = pymath.sqrt(2)
  nsssphe = math.nss_spherical_harmonics(order, 5000, lfg)
  for l in range(1, lmax+1):
    for m in range(-l, l+1):
      for it, t in enumerate(theta):
        for ip, p in enumerate(phi):
          Ylm = nsssphe.spherical_harmonic(l, abs(m), t, p)
          if m < 0:
            r = sqrt2 * ((-1) ** m) * Ylm.imag
          elif m == 0:
            assert Ylm.imag == 0.0
            r = Ylm.real
          else:
            r = sqrt2 * ((-1) ** m) * Ylm.real
          Intensity[ip, it] += params[counter] * r
      counter += 1
  Intensity = 1.0/Intensity #make scale factor, not inverse
  
  if Intensity.max() - Intensity.min() != 0.0:
    rel_Int = (Intensity - Intensity.min())/(Intensity.max() - Intensity.min())
  else:
    rel_Int = Intensity
  #print "max, min absorption factors are (%s,%s)" % (Intensity.max(),Intensity.min())
  plt.figure(figsize=(8,6))
  gs = gridspec.GridSpec(1, 1)
  ax = plt.subplot(gs[0, 0])
  im = ax.imshow(rel_Int.T, cmap='viridis', origin='lower')
  ax.set_yticks([0,  (STEPS-1)/4.0, (STEPS-1)/2.0, 3.0*(STEPS-1)/4.0, STEPS-1])
  ax.set_yticklabels([0, 45, 90, 135, 180])
  ax.set_ylabel('Phi (degrees)')
  ax.set_xticks([0.0, (2.0*STEPS-1)/4.0, (2.0*STEPS-1)/2.0, 3.0*(2.0*STEPS-1)/4.0, (2.0*STEPS-1)])
  ax.set_xticklabels([0, 90, 180, 270, 360])
  ax.set_xlabel('Theta (degrees)')
  ax.set_title('Scale factors for absorption correction (note: not inverse scales)')
  divider = make_axes_locatable(ax)
  cax1 = divider.append_axes("right", size="5%", pad=0.05)
  cbar = plt.colorbar(im, cax=cax1,ticks=[0, 0.25, 0.5, 0.75, 1])
  cbar.ax.set_yticklabels(['%.6f' % Intensity.min(), 
                           '%.6f' % ((Intensity.min()*0.75) + (Intensity.max()*0.25)), 
                           '%.6f' % ((Intensity.min()*0.5) + (Intensity.max()*0.5)),
                           '%.6f' % ((Intensity.min()*0.25) + (Intensity.max()*0.75)),
                           '%.6f' % Intensity.max()])
  if outputfile:
    plt.savefig(outputfile)
  else:
    plt.savefig('absorption_surface')

if __name__ == "__main__":
  datafile="/Users/whi10850/Documents/dials_scratch/jbe/scaling_code/test_data/x4_wide_integrated_scaled.pickle"
  data_man = load_data(filename = datafile)

  plot_data_decay(data_man)
  plot_data_absorption(data_man)
  plot_data_modulation(data_man)
