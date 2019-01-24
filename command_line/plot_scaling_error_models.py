from __future__ import absolute_import, division, print_function
import sys
import matplotlib
from libtbx.utils import Sorry
from dials.util import halraiser
from dials.util.options import OptionParser, flatten_reflections, flatten_experiments
from dials.array_family import flex
from dials.algorithms.scaling.Ih_table import IhTable
from dials.algorithms.scaling.error_model.error_model import \
  BasicErrorModel
from dials.algorithms.scaling.scaling_library import \
  create_scaling_model
from dials.algorithms.scaling.scaler_factory import create_scaler
from dials.algorithms.scaling.scaler import MultiScalerBase
from dials.util.multi_dataset_handling import assign_unique_identifiers, \
  parse_multiple_datasets
from libtbx import phil
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

phil_scope = phil.parse('''
  output {
    plot_out = "error_models.png"
      .type = str
      .help = "Option to set filename for error_model plot."
    combine_plots = False
      .type = bool
      .help = "Option to plot all datasets on one plot rather than individually"
    plot_labels = None
      .type = strings
  }
include scope dials.algorithms.scaling.scaling_options.phil_scope
''', process_includes=True)



def main(argv):
  '''the plotting script'''

  optionparser = OptionParser(usage=None, read_experiments=True,
    read_reflections=True, phil=phil_scope,
    check_format=False)
  params, _ = optionparser.parse_args(argv, show_diff_phil=False)


  if not params.input.reflections:
    optionparser.print_help()
    return

  diff_phil = optionparser.diff_phil.as_str()
  if diff_phil is not '':
    print('The following parameters have been modified:\n')
    print(diff_phil)

  reflections = flatten_reflections(params.input.reflections)
  experiments = flatten_experiments(params.input.experiments)

  colors = ['b', 'r', 'g', 'k'] * len(reflections)

  if len(experiments) != 1:
    print(('Checking for the existence of a reflection table {sep}'
      'containing multiple scaled datasets {sep}').format(sep='\n'))
    reflections = parse_multiple_datasets(reflections)
    print("Found %s reflection tables in total." % len(reflections))
    print("Found %s experiments in total." % len(experiments))

  if params.output.plot_labels:
    if len(params.output.plot_labels) != len(reflections):
      print("Length of plot labels not equal to number of datasets, exiting")
      sys.exit()

  experiments, reflections = assign_unique_identifiers(
      experiments, reflections)

  experiments = create_scaling_model(params, experiments, reflections)


  #for refl_table, exp in zip(reflections, experiments):
  scaler = create_scaler(params, experiments, reflections)
  if isinstance(scaler, MultiScalerBase):
    if params.output.combine_plots:
      plt.figure(figsize=(8, 8))
      gs = gridspec.GridSpec(2, 2)
    for i, single_scaler in enumerate(scaler.single_scalers):
      good_sel = ~(single_scaler.reflection_table.get_flags(
        single_scaler.reflection_table.flags.bad_for_scaling, all=False))
      single_scaler.Ih_table = IhTable([(single_scaler.reflection_table, good_sel)],
        scaler.space_group)
      Ih_table = single_scaler.Ih_table.blocked_data_list[0]
      error_model = BasicErrorModel(Ih_table)
      if not 'error_model_parameters' in single_scaler.experiments.scaling_model.configdict:
        raise Sorry("""No error model found in experiments file - likely cause is that
          no error model was optimised for this dataset.""")
      error_model.refined_parameters = single_scaler.experiments.scaling_model.configdict[
        'error_model_parameters']
      error_model.update_for_minimisation(error_model.refined_parameters)
      if params.output.plot_labels:
        label = params.output.plot_labels[i]
      else:
        label = None
      if params.output.combine_plots:
        normal_probability_plot(error_model, mode='add', filename='none', gs=gs,
          opacity=1.0/len(reflections), label=label, color=colors[i])
      else:
        normal_probability_plot(error_model, filename="error_models_"+str(i)+".png", label=label)
    if params.output.combine_plots:
      plt.tight_layout()
      filename = "combined_error_models.png"
      print("Saving plot to %s" % filename)
      plt.savefig(filename)
  else:
    good_sel = ~(scaler.reflection_table.get_flags(
      scaler.reflection_table.flags.bad_for_scaling, all=False))
    scaler.Ih_table = IhTable([(scaler.reflection_table, good_sel)], scaler.space_group)
    Ih_table = scaler.Ih_table.blocked_data_list[0]
    error_model = BasicErrorModel(Ih_table)
    if not 'error_model_parameters' in scaler.experiments.scaling_model.configdict:
      raise Sorry("""No error model found in experiments file - likely cause is that
        no error model was optimised for this dataset.""")
    error_model.refined_parameters = scaler.experiments.scaling_model.configdict[
      'error_model_parameters']
    error_model.update_for_minimisation(error_model.refined_parameters)
    if params.output.plot_labels:
      label = params.output.plot_labels[0]
    else:
      label = None
    normal_probability_plot(error_model, filename=params.output.plot_out, label=label)

def normal_probability_plot(error_model, filename, mode='full', gs='none',
  opacity=1.0, label=None, color='b'):
  from scitbx.math import distributions
  import numpy as np
  norm = distributions.normal_distribution()

  n = len(error_model.delta_hl)
  if n <= 10:
    a = 3/8
  else:
    a = 0.5

  sorted_data = flex.sorted(error_model.delta_hl)
  rankits = flex.double([norm.quantile((i+1-a)/(n+1-(2*a))) for i in xrange(n)])
  if mode == 'full':
    plt.figure(figsize=(8, 8))
    gs = gridspec.GridSpec(2, 2)
  ax1 = plt.subplot(gs[0, 0])
  ax2 = plt.subplot(gs[0, 1])
  ax3 = plt.subplot(gs[1, 0])
  ax4 = plt.subplot(gs[1, 1])
  x = np.linspace(-4,4,100)
  ax1.scatter(rankits, sorted_data, alpha=opacity, label=label, color=color, s=2)
  ax1.set_title("Normal probability plot with error model applied")
  ax1.set_xlabel("Order statistic medians, m")
  ax1.set_ylabel("Ordered responses, z")
  ax1.plot(x, x, color='k')
  ax2.set_xlabel("Normalised deviation (z)")
  ax2.set_ylabel("Counts")
  ax2.set_title("z histogram with error model applied")
  ax2.hist(sorted_data, bins=100, alpha=opacity, label=label, color=color)
  ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  error_model.refined_parameters = [1.0, 0.0]
  error_model.update_for_minimisation(error_model.refined_parameters)
  sorted_data = flex.sorted(error_model.delta_hl)

  ax3.scatter(rankits, sorted_data, alpha=opacity, label=label, color=color, s=2)
  ax3.plot(x, x, color='k')
  ax3.set_title("Normal probability plot without error model applied")
  ax3.set_xlabel("Order statistic medians, m")
  ax3.set_ylabel("Ordered responses, z")
  ax4.set_xlabel("Normalised deviation (z) ")
  ax4.set_ylabel("Counts")
  ax4.set_title("z histogram without error model applied")
  ax4.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  ax4.hist(sorted_data, bins=100, alpha=opacity, label=label, color=color)
  if label:
    ax1.legend(loc=2)
    ax2.legend(loc=2)
    ax3.legend(loc=2)
    ax4.legend(loc=2)
  if mode == 'full':
    plt.tight_layout()
    print("Saving plot to %s" % filename)
    plt.savefig(filename)


if __name__ == "__main__":
  try:
    sys.exit(main(sys.argv[1:]))
  except Exception as e:
    halraiser(e)
