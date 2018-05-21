from __future__ import absolute_import, division, print_function
import sys
import matplotlib
from libtbx.utils import Sorry
from dials.util import halraiser
from dials.util.options import OptionParser, flatten_reflections, flatten_experiments
from dials.array_family import flex
from dials.algorithms.scaling.Ih_table import SingleIhTable
from dials.algorithms.scaling.error_model.error_model import \
  BasicErrorModel
from dials.algorithms.scaling.model.scaling_model_factory import \
  create_scaling_model
from dials.algorithms.scaling.scaler_factory import create_scaler
from dials.algorithms.scaling.scaler import MultiScalerBase
from dials.algorithms.scaling.scaling_utilities import parse_multiple_datasets
from libtbx import phil
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

phil_scope = phil.parse('''
  output {
    plot_out = "error_models.png"
      .type = str
      .help = "Option to set filename for error_model plot."
  }
include scope dials.algorithms.scaling.scaling_options.phil_scope
''', process_includes=True)

def main(argv):
  '''the plotting script'''

  optionparser = OptionParser(usage=None, read_experiments=True,
    read_reflections=True, read_datablocks=False, phil=phil_scope,
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

  if len(experiments) != 1:
    print(('Checking for the existence of a reflection table {sep}'
      'containing multiple scaled datasets {sep}').format(sep='\n'))
    reflections = parse_multiple_datasets(reflections)
    print("Found %s reflection tables in total." % len(reflections))
    print("Found %s experiments in total." % len(experiments))

  experiments = create_scaling_model(params, experiments, reflections)


  #for refl_table, exp in zip(reflections, experiments):
  scaler = create_scaler(params, experiments, reflections)
  if isinstance(scaler, MultiScalerBase):
    for i, single_scaler in enumerate(scaler.single_scalers):
      single_scaler.expand_scales_to_all_reflections()
      weighting = single_scaler._update_weights_for_scaling(single_scaler.reflection_table,
        single_scaler.params, weights_filter=False, error_model_params=None)
      single_scaler.Ih_table = SingleIhTable(single_scaler.reflection_table, scaler.space_group,
        weighting.weights)
      Ih_table = single_scaler.Ih_table
      error_model = BasicErrorModel(Ih_table)
      if not 'error_model_parameters' in single_scaler.experiments.scaling_model.configdict:
        raise Sorry("""No error model found in experiments file - likely cause is that
          no error model was optimised for this dataset.""")
      error_model.refined_parameters = single_scaler.experiments.scaling_model.configdict[
        'error_model_parameters']
      error_model.update_for_minimisation(error_model.refined_parameters)
      normal_probability_plot(error_model, filename="error_models_"+str(i)+".png")
  else:
    scaler.expand_scales_to_all_reflections()
    weighting = scaler._update_weights_for_scaling(scaler.reflection_table, scaler.params,
      weights_filter=False, error_model_params=None)
    scaler.Ih_table = SingleIhTable(scaler.reflection_table, scaler.space_group, weighting.weights)
    Ih_table = scaler.Ih_table
    error_model = BasicErrorModel(Ih_table)
    if not 'error_model_parameters' in scaler.experiments.scaling_model.configdict:
      raise Sorry("""No error model found in experiments file - likely cause is that
        no error model was optimised for this dataset.""")
    error_model.refined_parameters = scaler.experiments.scaling_model.configdict[
      'error_model_parameters']
    error_model.update_for_minimisation(error_model.refined_parameters)
    normal_probability_plot(error_model, filename=params.output.plot_out)

def normal_probability_plot(error_model, filename):
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
  from matplotlib import pyplot as plt
  plt.figure(figsize=(8, 8))
  gs = gridspec.GridSpec(2, 2)
  ax1 = plt.subplot(gs[0, 0])
  ax2 = plt.subplot(gs[0, 1])
  ax3 = plt.subplot(gs[1, 0])
  ax4 = plt.subplot(gs[1, 1])
  x = np.linspace(-4,4,100)
  ax1.scatter(sorted_data, rankits)
  ax1.set_title("Normal probability plot, errors refined.")
  ax1.set_xlabel("Distribution of data")
  ax1.set_ylabel("Theoretical normal distribution")
  ax1.plot(x, x)
  ax2.set_xlabel("Normalised deviation from the mean")
  ax2.set_ylabel("Count")
  ax2.set_title("Data")
  ax2.hist(sorted_data, bins=100)

  error_model.refined_parameters = [1.0, 0.0]
  error_model.update_for_minimisation(error_model.refined_parameters)
  sorted_data = flex.sorted(error_model.delta_hl)

  ax3.scatter(sorted_data, rankits)
  ax3.plot(x, x)
  ax3.set_title("Normal probability plot, errors unrefined.")
  ax3.set_xlabel("Distribution of data")
  ax3.set_ylabel("Theoretical normal distribution")
  ax4.set_xlabel("Normalised deviation from the mean")
  ax4.set_ylabel("Count")
  ax4.set_title("Data")
  ax4.hist(sorted_data, bins=100)
  plt.tight_layout()
  print("Saving plot to %s" % filename)
  plt.savefig(filename)


if __name__ == "__main__":
  try:
    sys.exit(main(sys.argv[1:]))
  except Exception as e:
    halraiser(e)
