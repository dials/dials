"""
Plot the scaling models from a scaled_experiments.json and scaled.pickle
"""
from __future__ import absolute_import, division, print_function
import sys
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from dials.array_family import flex
from dials.util import halraiser
from dials.util.options import OptionParser, flatten_reflections, flatten_experiments
from dials.algorithms.scaling.model import model as Model
from dials.algorithms.scaling.scaling_library import create_scaling_model
from dials.algorithms.scaling.scaling_utilities import parse_multiple_datasets
from libtbx import phil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


phil_scope = phil.parse('''
  debug = False
    .type = bool
    .help = "Output additional debugging information"
  output {
    scales_out = "scale_model"
      .type = str
      .help = "Option to set filename for output plot."
    absorption_out = "absorption_surface"
      .type = str
      .help = "Option to set filename for output absorption plot."
    decay_out = "decay_correction"
      .type = str
      .help = "Option to set filename for 2D decay plot
               (only relevant for array-based parameterisation)."
    abscorplot = "absorption_correction"
      .type = str
      .help = "Option to set filename for absorption plot
               (only relevant for array-based parameterisation)."
    modcorplot = "modulation_correction"
      .type = str
      .help = "Option to set filename for modulation plot
               (only relevant for array-based parameterisation)."
    experiment_range = ()
      .type = ints
      .help = "Option to specify a subset of experiments to plot."
    with_errors = True
      .type = bool
      .help = "Option to turn off plotting of parameters and errors."
    limit_range_to_obs = False
      .type = bool
      .help = "Option to only plot parameters and errors within experimental range."
  }
  include scope dials.algorithms.scaling.scaling_options.phil_scope
''', process_includes=True)

def plot_scaling_models(argv):
  '''the plotting script'''

  optionparser = OptionParser(usage=None, read_experiments=True,
    read_reflections=True, read_datablocks=False, phil=phil_scope,
    check_format=False)
  params, _ = optionparser.parse_args(argv, show_diff_phil=False)


  if not params.input.experiments or not params.input.reflections:
    optionparser.print_help()
    return

  diff_phil = optionparser.diff_phil.as_str()
  if diff_phil is not '':
    print('The following parameters have been modified:\n')
    print(diff_phil)

  reflections = flatten_reflections(params.input.reflections)
  experiments = flatten_experiments(params.input.experiments)

  if len(experiments) != 1:
    print(('Checking for the existence of a reflection table containing {sep}'
    'multiple scaled datasets. {sep}').format(sep='\n'))
    reflections, ids = parse_multiple_datasets(reflections)
    print("Found %s experiments in total." % len(experiments))
    print('\n'+'*'*40)

  experiments = create_scaling_model(params, experiments, reflections)

  print('\nPlotting graphs of scale factors. \n')
  if len(experiments) == 1:
    experiment = experiments[0]
    reflection = reflections[0]
    if isinstance(experiment.scaling_model, Model.PhysicalScalingModel):
      if ('scale' in experiment.scaling_model.configdict['corrections'] or
        'decay' in experiment.scaling_model.configdict['corrections']):
        plot_smooth_scales(params, experiment, reflection,
          outputfile=str(params.output.scales_out)+'.png')
      if 'absorption' in experiment.scaling_model.configdict['corrections']:
        plot_absorption_surface(experiment,
          outputfile=str(params.output.absorption_out)+'.png')
    elif isinstance(experiment.scaling_model, Model.ArrayScalingModel):
      if 'decay' in experiment.scaling_model.configdict['corrections']:
        plot_2D_decay_correction(experiment, reflection)
      if 'absorption' in experiment.scaling_model.configdict['corrections']:
        plot_3D_absorption_correction(experiment, reflection)
      if 'modulation' in experiment.scaling_model.configdict['corrections']:
        plot_2D_modulation_correction(experiment, reflection)
  else:
    if params.output.experiment_range:
      for j in params.output.experiment_range:
        experiment = experiments[j]
        reflection = reflections[j]
        plot_multi(params, experiment, reflection, j)
    else:
      for j, (experiment, reflection) in enumerate(zip(experiments, reflections)):
        plot_multi(params, experiment, reflection, j)
  print("\nFinished plotting graphs of scale factors. \n")


def plot_multi(params, experiment, reflection, j):
  '''subscript to plot a single instance of a multi-dataset file'''
  if isinstance(experiment.scaling_model, Model.PhysicalScalingModel):
    if ('scale' in experiment.scaling_model.configdict['corrections'] or
      'decay' in experiment.scaling_model.configdict['corrections']):
      plot_smooth_scales(params, experiment, reflection,
        outputfile=str(params.output.scales_out)+'_'+str(j+1)+'.png')
    if 'absorption' in experiment.scaling_model.configdict['corrections']:
      plot_absorption_surface(experiment,
        outputfile=str(params.output.absorption_out)+'_'+str(j+1)+'.png')
  elif isinstance(experiment.scaling_model, Model.ArrayScalingModel):
    if 'decay' in experiment.scaling_model.configdict['corrections']:
      plot_2D_decay_correction(experiment, reflection,
        outputfile=str(params.output.decay_out)+'_'+str(j+1)+'.png')
    if 'absorption' in experiment.scaling_model.configdict['corrections']:
      plot_3D_absorption_correction(experiment, reflection,
        outputfile=str(params.output.abscorplot)+'_'+str(j+1)+'.png')
    if 'modulation' in experiment.scaling_model.configdict['corrections']:
      plot_2D_modulation_correction(experiment, reflection,
        outputfile=str(params.output.modcorplot)+'_'+str(j+1)+'.png')

def plot_smooth_scales(params, experiments, reflections, outputfile=None):
  """Plot the scale and decay terms of a physical scaling model."""
  plt.figure(figsize=(8.5, 5))
  legends = []
  ax1 = None
  reflections = reflections.select(~reflections.get_flags(
    reflections.flags.user_excluded_in_scaling))
  reflections = reflections.select(reflections.get_flags(
    reflections.flags.integrated))
  if 'scale' in experiments.scaling_model.configdict['corrections']:
    reflections['norm_rot_angle'] = (reflections['xyzobs.px.value'].parts()[2]
      * experiments.scaling_model.configdict['s_norm_fac'])
    reflections['norm_rot_angle'] = (reflections['norm_rot_angle']
      - min(reflections['norm_rot_angle']))
    scale_rot_int = experiments.scaling_model.configdict['scale_rot_interval']
    int_rel_max = int(max(reflections['norm_rot_angle'])) + 1
    int_rel_min = (int(min(reflections['norm_rot_angle'])))
    rel_values = flex.double(np.linspace(0, int_rel_max-int_rel_min,
      ((int_rel_max-int_rel_min)/0.1)+1, endpoint=True))
    rel_values[-1] = rel_values[-1] - 0.0001
    rt = flex.reflection_table()
    rt['norm_rot_angle'] = rel_values
    #rel_values = rel_values - rel_values[0]
    scale_SF = experiments.scaling_model.components['scale']
    scale_SF.update_reflection_data(rt)
    scale_SF.calculate_scales()
    smoother_phis = [i * scale_rot_int for i in scale_SF.smoother.positions()]
    ax1 = plt.subplot(2, 1, 1)
    #plt.title('Smooth scale factors')
    ax1.plot(rel_values*scale_rot_int, scale_SF.inverse_scales[0],
      label='smootly varying \ninverse scale factor')
    if params.output.with_errors:
      if params.output.limit_range_to_obs:
        ax1.errorbar(smoother_phis[1:-1], scale_SF.parameters[1:-1],
          yerr=scale_SF.parameter_esds[1:-1], fmt='o', color='k')
      else:
        ax1.errorbar(smoother_phis, scale_SF.parameters,
          yerr=scale_SF.parameter_esds, fmt='o', color='k',
          label='model parameters')
    ax1.set_ylabel('Inverse scale factor', fontsize=12)
    ax1.set_xlabel('Rotation angle ('+r'$^{\circ}$'+')', fontsize=12)
    ax1.tick_params(axis='both', which='major', labelsize=12)
    leg1 = ax1.legend(bbox_to_anchor=(1.02, 1), loc=2, fontsize=12)
    legends.append(leg1)

  if 'decay' in experiments.scaling_model.configdict['corrections']:
    reflections['norm_time_values'] = (reflections['xyzobs.px.value'].parts()[2]
      * experiments.scaling_model.configdict['d_norm_fac'])
    reflections['norm_time_values'] = (reflections['norm_time_values']
      - min(reflections['norm_time_values']))
    decay_rot_int = experiments.scaling_model.configdict['decay_rot_interval']
    int_rel_max = int(max(reflections['norm_time_values'])) + 1
    int_rel_min = (int(min(reflections['norm_time_values'])))
    rel_values = flex.double(np.linspace(0, int_rel_max-int_rel_min,
      ((int_rel_max-int_rel_min)/0.1)+1, endpoint=True))
    rel_values[-1] = rel_values[-1] - 0.0001
    #rel_values = rel_values# - rel_values[0]
    decay_SF = experiments.scaling_model.components['decay']
    rt = flex.reflection_table()
    rt['norm_time_values'] = rel_values
    rt['d'] = flex.double(rel_values.size(), 1.0)
    decay_SF.update_reflection_data(rt)#normalised_values=rel_values,
    #  dvalues=flex.double(rel_values.size(), 1.0))
    decay_SF.calculate_scales()
    smoother_phis = [i * decay_rot_int for i in decay_SF._smoother.positions()]
    if ax1:
      ax2 = plt.subplot(2, 1, 2, sharex=ax1)
    else:
      ax2 = plt.subplot(2, 1, 2)
    ax2.set_ylabel('Relative B factor (' + r'$\AA^{2}$'+')', fontsize=12)
    ax2.set_xlabel('Rotation angle (' + r'$^{\circ}$'+')', fontsize=12)
    ax2.plot(rel_values * decay_rot_int, np.log(decay_SF.inverse_scales[0])*2.0,
      label='smootly varying \nB-factor') #convert scales to B values
    if params.output.with_errors:
      if params.output.limit_range_to_obs:
        ax2.errorbar(smoother_phis[1:-1], decay_SF.parameters[1:-1],
          yerr=decay_SF.parameter_esds[1:-1], fmt='o', color='k')
      else:
        ax2.errorbar(smoother_phis, decay_SF.parameters,
          yerr=decay_SF.parameter_esds, fmt='o', color='k', label='model parameters')
    leg2 = ax2.legend(bbox_to_anchor=(1.02, 1), loc=2, fontsize=12)
    ax2.tick_params(axis='both', which='major', labelsize=12)
    if ax1:
      plt.setp(ax1.get_xticklabels(), visible=False)
      ax1.set_xlabel('')
    legends.append(leg2)
  plt.tight_layout()
  if outputfile:
    plt.savefig(outputfile, bbox_extra_artists=(legends), bbox_inches='tight')
    print("Saving plot to %s" % outputfile)
  else:
    plt.savefig('smooth_scale_factors.png')
    print("Saving plot to smooth_scale_factors.png")

def plot_absorption_surface(experiment, outputfile=None):
  """Plot a spherical harmonic absorption surface."""
  params = experiment.scaling_model.components['absorption'].parameters

  from scitbx import math
  import math as pymath
  order = int(-1.0 + ((1.0 + len(params))**0.5))
  lfg = math.log_factorial_generator(2 * order + 1)
  STEPS = 50
  phi = np.linspace(0, 2 * np.pi, 2*STEPS)
  theta = np.linspace(0, np.pi, STEPS)
  THETA, _ = np.meshgrid(theta, phi)
  lmax = int(-1.0 + ((1.0 + len(params))**0.5))
  Intensity = np.ones(THETA.shape)
  counter = 0
  sqrt2 = pymath.sqrt(2)
  nsssphe = math.nss_spherical_harmonics(order, 50000, lfg)
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
  #Intensity = 1.0/Intensity #make scale factor, not inverse

  if Intensity.max() - Intensity.min() != 0.0:
    rel_Int = (Intensity - Intensity.min())/(Intensity.max() - Intensity.min())
  else:
    rel_Int = Intensity
  #print "max, min absorption factors are (%s,%s)" % (Intensity.max(),Intensity.min())
  plt.figure(figsize=(8, 5))
  gs = gridspec.GridSpec(1, 1)
  ax = plt.subplot(gs[0, 0])
  im = ax.imshow(rel_Int.T, cmap='viridis', origin='lower')
  ax.set_yticks([0, (STEPS-1)/4.0, (STEPS-1)/2.0, 3.0*(STEPS-1)/4.0, STEPS-1])
  ax.set_yticklabels([0, 45, 90, 135, 180], fontsize=12)
  ax.set_ylabel(r'$\phi$'+' ('+r'$^{\circ}$'+')', fontsize=14)
  ax.set_xticks([0.0, (2.0*STEPS-1)/4.0, (2.0*STEPS-1)/2.0, 3.0*(2.0*STEPS-1)/4.0, (2.0*STEPS-1)])
  ax.set_xticklabels([0, 90, 180, 270, 360], fontsize=12)
  ax.set_xlabel(r'$\theta$'+' ('+r'$^{\circ}$'+')', fontsize=14)
  ax.set_title('Absorption correction surface - inverse scale factors', fontsize=12)
  divider = make_axes_locatable(ax)
  cax1 = divider.append_axes("right", size="5%", pad=0.05)
  cbar = plt.colorbar(im, cax=cax1, ticks=[0, 0.25, 0.5, 0.75, 1])
  cbar.ax.set_yticklabels(['%.3f' % Intensity.min(),
                           '%.3f' % ((Intensity.min()*0.75) + (Intensity.max()*0.25)),
                           '%.3f' % ((Intensity.min()*0.5) + (Intensity.max()*0.5)),
                           '%.3f' % ((Intensity.min()*0.25) + (Intensity.max()*0.75)),
                           '%.3f' % Intensity.max()], fontsize=12)
  #ax2 = plt.subplot(gs[1, 0])
  #data = [list(params)[0:5],list(params)[5:10],list(params)[10:15],list(params)[15:20]]
  #errors = [list(experiment.scaling_model.components['absorption'].parameter_esds)]
  #ax2.table(cellText=data, loc='top')
  #ax2.axis("off")
  plt.tight_layout()
  if outputfile:
    plt.savefig(outputfile)
    print("Saving plot to %s" % outputfile)
  else:
    plt.savefig('absorption_surface.png')
    print("Saving plot to absorption_surface.png")

def plot_2D_decay_correction(experiment, reflections, outputfile=None):
  '''plotting of decay vs time correction for array-based parameterisation'''
  '''first extract the model and data'''
  reflections = reflections.select(reflections['d'] > 0.0)
  reflections = reflections.select(~reflections.get_flags(
    reflections.flags.user_excluded_in_scaling))
  reflections = reflections.select(reflections.get_flags(
    reflections.flags.integrated))
  configdict = experiment.scaling_model.configdict
  reflections['normalised_res_values'] = (((1.0 / (reflections['d']**2))
      - configdict['resmin']) / configdict['res_bin_width'])
  reflections['norm_time_values'] = (reflections['xyzobs.px.value'].parts()[2]
      * configdict['time_norm_fac'])
  reflections['norm_time_values'] = (reflections['norm_time_values']
      - min(reflections['norm_time_values']))

  time_rot_int = configdict['time_rot_interval']
  int_rel_max = int(max(reflections['norm_time_values'])) + 1
  int_rel_min = (int(min(reflections['norm_time_values'])))

  rel_values = flex.double(np.linspace(0, int_rel_max - int_rel_min,
    ((int_rel_max-int_rel_min)/0.1)+1, endpoint=True))
  rel_values[-1] = rel_values[-1] - 0.0001
  rel_values_2 = rel_values - rel_values[0]
  #x_axis_vals = rel_values_2

  '''create a grid of x and y points and use these to generate scale factors'''
  max_res = int(max(reflections['normalised_res_values'])) + 1
  #max_time = int(max(reflections['norm_time_values'])) + 1
  rel_values_1 = np.arange(0, max_res+0.1, 0.1)
  #rel_values_2 = np.arange(0, max_time+0.1, 0.1)
  (n1, n2) = (len(rel_values_1), len(rel_values_2))
  rel_values_1 = flex.double(np.tile(rel_values_1, n2))
  rel_values_2 = flex.double(np.repeat(rel_values_2, n1))

  rt = flex.reflection_table()
  rt['norm_time_values'] = rel_values_2
  rt['normalised_res_values'] = rel_values_1

  decay_factor = experiment.scaling_model.components['decay']
  decay_factor.update_reflection_data(rt)
  decay_factor.calculate_scales()
  scales = decay_factor.inverse_scales
  scalefactor_2D = np.reshape(list(scales), (n2, n1)).T

  '''generate a plot'''
  plt.figure(figsize=(8, 4))
  gs = gridspec.GridSpec(1, 2)
  ax1 = plt.subplot(gs[0, 0])
  ax2 = plt.subplot(gs[0, 1])

  resmax = (1.0 / (min(reflections['d'])**2))
  resmin = (1.0 / (max(reflections['d'])**2))
  resbin_boundaries = np.arange(resmin, resmax, 2*(resmax - resmin)/(
    experiment.scaling_model.configdict['n_res_param']-1))
  dbin_boundaries = ['%.3f' % x for x in 1.0/(resbin_boundaries**0.5)]
  im = ax1.imshow(scalefactor_2D, cmap='viridis', origin='upper', aspect='auto')
  divider = make_axes_locatable(ax1)
  cax1 = divider.append_axes("right", size="5%", pad=0.05)
  _ = plt.colorbar(im, cax=cax1)
  ax1.set_ylabel('d ('+r'$\AA$'+')')
  ax1.set_xlabel('Rotation angle ('+r'$^{\circ}$'+')')
  ax1.set_yticks(np.arange(0, (max_res * 10)+0.01, 20))
  ax1.set_yticklabels(dbin_boundaries)
  #ax1.set_xticks(x_axis_vals)
  #ax1.set_xticklabels(x_axis_vals*time_rot_int)
  ax1.set_xticks(np.arange(0, ((int_rel_max - int_rel_min) * 10)+0.01, 10))
  xlabels = np.arange(0, (time_rot_int*(int_rel_max-int_rel_min))+0.001, time_rot_int)
  xlabels = np.around(xlabels, 1)
  ax1.set_xticklabels(xlabels)

  ax1.set_title('Decay correction (inverse scale factors)\n', fontsize=10)

  '''recalculate scales for plotting distribution in dataset'''
  decay_factor.update_reflection_data(rt)
  decay_factor.calculate_scales()

  ax2.hist(list(decay_factor.inverse_scales), 40, log=False)
  ax2.set_xlabel('inverse scale factor')
  ax2.set_ylabel('Counts')
  ax2.set_title('Distribution of dataset reflection corrections\n', fontsize=10)

  plt.tight_layout()
  if outputfile:
    plt.savefig(outputfile)
    print("Saving plot to %s" % outputfile)
  else:
    plt.savefig('decay_correction.png')
    print("Saving plot to decay_correction.png")

def plot_2D_modulation_correction(experiment, reflections, outputfile=None):
  """Plot the detector modulation correction for an array model."""
  configdict = experiment.scaling_model.configdict
  nxdet = ((reflections['xyzobs.px.value'].parts()[0] - configdict['xmin']) /
    configdict['x_det_bin_width'])
  nydet = ((reflections['xyzobs.px.value'].parts()[1] - configdict['ymin']) /
    configdict['y_det_bin_width'])


  #nxbins = data_man.g_modulation.n1_parameters
  #nybins = data_man.g_modulation.n2_parameters
  max_x = int(max(nxdet)) + 1
  max_y = int(max(nydet)) + 1
  '''create a grid of x and y points and use these to generate scale factors'''
  rel_values_1 = np.arange(0, max_x+0.2, 0.2)
  rel_values_2 = np.arange(0, max_y+0.2, 0.2)
  (n1, n2) = (len(rel_values_1), len(rel_values_2))
  rel_values_1 = np.tile(rel_values_1, n2)
  rel_values_2 = np.repeat(rel_values_2, n1)

  rt = flex.reflection_table()
  rt['normalised_y_det_values'] = rel_values_2
  rt['normalised_x_det_values'] = rel_values_1

  modulation_factor = experiment.scaling_model.components['decay']
  modulation_factor.update_reflection_data(rt)
  modulation_factor.calculate_scales()
  scales = modulation_factor.inverse_scales
  scalefactor_2D = np.reshape(list(scales), (n2, n1))

  plt.figure(figsize=(8, 4))
  gs = gridspec.GridSpec(1, 2)
  ax1 = plt.subplot(gs[0, 0])
  ax2 = plt.subplot(gs[0, 1])

  im = ax1.imshow(scalefactor_2D, cmap='viridis', origin='lower')
  divider = make_axes_locatable(ax1)
  cax1 = divider.append_axes("right", size="5%", pad=0.05)
  _ = plt.colorbar(im, cax=cax1)
  ax1.set_yticks(np.arange(0, (max_x * 5)+0.01, 10))
  ax1.set_yticklabels(np.arange(0, max_x+0.01, 2))
  ax1.set_xticks(np.arange(0, (max_y * 5)+0.01, 10))
  ax1.set_xticklabels(np.arange(0, max_y+0.01, 2))
  ax1.set_ylabel('Detector normalised y value')
  ax1.set_xlabel('Detector normalised x value')
  ax1.set_title('Map of detector correction (inverse scale factors)', fontsize=10)

  '''recalculate scales for plotting distribution in dataset'''
  modulation_factor.update_reflection_data(rt)
  modulation_factor.calculate_scales()

  ax2.hist(list(modulation_factor.inverse_scales), 40, log=False)
  ax2.set_xlabel('Detector correction inverse scale factor')
  ax2.set_ylabel('Counts')
  ax2.set_title('Distribution of dataset reflection corrections', fontsize=10)

  plt.tight_layout()
  if outputfile:
    plt.savefig(outputfile)
    print("Saving plot to %s" % outputfile)
  else:
    plt.savefig('modulation_correction.png')
    print("Saving plot to modulation_correction.png")

def plot_3D_absorption_correction(experiment, reflections, outputfile=None):
  '''produces a plot of the absorption correction model parameters
  for array-based parameterisation'''
  '''first extract the model and data'''
  configdict = experiment.scaling_model.configdict
  n_time_bins = configdict['n_time_param']
  n_abs_bins = configdict['n_x_param'] * configdict['n_y_param']
  (xvalues, yvalues, zvalues) = reflections['xyzobs.px.value'].parts()
  '''calculate normalised x, y, time values'''
  nt = zvalues * experiment.scaling_model.configdict['time_norm_fac']
  nax = (xvalues - configdict['xmin']) / configdict['x_bin_width']
  nay = (yvalues - configdict['ymin']) / configdict['y_bin_width']

  rt = flex.reflection_table()
  rt['normalised_x_abs_values'] = nax
  rt['normalised_y_abs_values'] = nay
  rt['norm_time_values'] = nt

  absorption_factor = experiment.scaling_model.components['absorption']
  absorption_factor.update_reflection_data(rt)
  absorption_factor.calculate_scales()
  parameters_2D = np.reshape(list(absorption_factor.parameters),
    (n_time_bins, n_abs_bins)).T

  '''generate a plot'''
  plt.figure(figsize=(8, 4))
  gs = gridspec.GridSpec(1, 2)
  ax1 = plt.subplot(gs[0, 0])
  ax2 = plt.subplot(gs[0, 1])

  im = ax1.imshow(parameters_2D, cmap='viridis', origin='lower')
  divider = make_axes_locatable(ax1)
  cax1 = divider.append_axes("right", size="5%", pad=0.05)
  _ = plt.colorbar(im, cax=cax1)
  ax1.set_ylabel('Detector position (n x n grid)')
  ax1.set_xlabel('Normalised time value')
  ax1.set_title('Model parameters for absorption correction', fontsize=10)

  ax2.hist(list(absorption_factor.inverse_scales), 40, log=False)
  ax2.set_xlabel('Absorption correction inverse scale factor')
  ax2.set_ylabel('Counts')
  ax2.set_title('Distribution of dataset reflection corrections', fontsize=10)

  plt.tight_layout()
  if outputfile:
    plt.savefig(outputfile)
    print("Saving plot to %s" % outputfile)
  else:
    plt.savefig('absorption_correction.png')
    print("Saving plot to absorption_correction.png")


if __name__ == "__main__":
  try:
    sys.exit(plot_scaling_models(sys.argv[1:]))
  except Exception as e:
    halraiser(e)
