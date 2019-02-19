"""
This program runs several rounds of scaling and filtering for multi-crystal
analysis, also producing summary plots.

This program runs consecutive cycles of dials.scale and dials.compute_delta_cchalf,
until one of several termination criteria are met - these are max_cycles,
max_percent_removed and min_completeness. Filtering is based on the calculation
of delta-cc-half values, that is the cc-half of the dataset when a subset of
images are removed. Any images with a highly negative cc-half are contributing
poorly to the overall merging statistics and are removed and the dataset is
then rescaled before further filtering analysis.

This program only acts to run the underlying programs and collect results in a
convenient manner - equivalent results should be obtained by running the
individual programs. All standard command-line options for these programs can
be given, alongside the termination criteria and output options for this program.
"""

from __future__ import absolute_import, division, print_function

import logging
import json
from collections import OrderedDict
from math import floor
import numpy as np
import matplotlib
import libtbx.phil
from cctbx.sgtbx import uctbx
import dials.util
from dials.command_line.scale import Script as ScalingScript
from dials.command_line.compute_delta_cchalf import Script as FilterScript
from dials.util.exclude_images import get_valid_image_ranges
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

# Define a logger
logger = logging.getLogger('dials')

# Define the master PHIL scope for this program
phil_scope = libtbx.phil.parse('''
filtering {
  max_cycles = 6
    .type = int(value_min=1)
  max_percent_removed = 10
    .type = float
  min_completeness = None
    .type = float(value_min=0, value_max=100)
    .help = "Desired minimum completeness, as a percentage (0 - 100)."
}
scale {
  include scope dials.command_line.scale.phil_scope
}
delta_cc_half {
  include scope dials.command_line.compute_delta_cchalf.phil_scope
}
output {
    log = dials.scale.log
      .type = str
      .help = "The log filename"
    debug.log = dials.scale.debug.log
      .type = str
      .help = "The debug log filename"
    experiments = "scaled_experiments.json"
      .type = str
      .help = "Option to set filepath for output json."
    reflections = "scaled.pickle"
      .type = str
      .help = "Option to set filepath for output pickle file of scaled
               intensities."
    analysis_results = "analysis_results.json"
      .type = str
      .help = "Option to set filepath for output json of analysis results."
    plots {
      histogram = cc_half_histograms.png
        .type = str
      merging_stats = merging_stats.png
        .type = str
      image_ranges = reduced_image_ranges.png
        .type = str
    }
}
''', process_includes=True)

class AnalysisResults(object):
  """Class to store results from scaling and filtering."""

  def __init__(self):
    self.termination_reason = None
    self.cycle_results = []
    self.initial_n_reflections = None
    self.initial_expids_and_image_ranges = None
    self.final_stats = None

  def add_cycle_results(self, results_dict):
    """Add the results dict from a scale and filter cycle."""
    merging_stats_dict = self._parse_merging_stats(results_dict['merging_stats'])
    results_dict['merging_stats'] = merging_stats_dict
    self.cycle_results.append(results_dict)

  @staticmethod
  def _parse_merging_stats(merging_stats_obj):
    merging_stats = {}
    overall = merging_stats_obj.overall
    merging_stats['ccs'] = [b.cc_one_half for b in merging_stats_obj.bins]
    merging_stats['rmerge'] = [b.r_merge for b in merging_stats_obj.bins]
    merging_stats['rpim'] = [b.r_pim for b in merging_stats_obj.bins]
    merging_stats['d_min'] = [b.d_min for b in merging_stats_obj.bins]
    merging_stats['overall'] = OrderedDict([('cc_one_half', overall.cc_one_half),
      ('r_merge', overall.r_merge), ('r_pim', overall.r_pim),
      ('i_over_sigma_mean', overall.i_over_sigma_mean),
      ('completeness', 100 * overall.completeness), ('n_obs', overall.n_obs)])
    return merging_stats

  def get_cycle_results(self):
    """Get the results from all cycles."""
    return self.cycle_results

  def get_last_cycle_results(self):
    """Get the results from the latest recorded cycle."""
    return self.cycle_results[-1]

  def add_final_stats(self, final_stats):
    """Add additional final merging stats from final rescale."""
    self.final_stats = self._parse_merging_stats(final_stats)

  def get_merging_stats(self):
    """Get all merging stats, including additional final stats if present."""
    stats = [res['merging_stats'] for res in self.cycle_results]
    if self.final_stats:
      stats += [self.final_stats]
    return stats

  def finish(self, termination_reason):
    """Set the termination reason/"""
    assert termination_reason in ['no_more_removed', 'max_cycles',
      'max_percent_removed', 'below_completeness_limit']
    self.termination_reason = termination_reason

  def to_dict(self):
    """Return the stored data as a dictionary."""
    return {'termination_reason' : self.termination_reason,
      'initial_n_reflections' : self.initial_n_reflections,
      'initial_expids_and_image_ranges' : self.initial_expids_and_image_ranges,
      'cycle_results' : OrderedDict([(i+1, val) for i, val in enumerate(self.cycle_results)]),
      'final_stats' : self.final_stats}

  @staticmethod
  def from_dict(dictionary):
    """Configure the class from its dictionary form."""
    results = AnalysisResults()
    results.termination_reason = dictionary['termination_reason']
    results.initial_expids_and_image_ranges = dictionary['initial_expids_and_image_ranges']
    results.cycle_results = [res for res in dictionary['cycle_results']]
    results.initial_n_reflections = dictionary['initial_n_reflections']
    results.final_stats = dictionary['final_stats']
    return results

class ScaleAndFilter(object):

  """Class to encapsulate a scaling and filtering algorithm."""

  def __init__(self, scaling_script, filtering_script):
    """Provide script classes that do scaling and filtering"""
    self.scaling_script = scaling_script
    self.filtering_script = filtering_script

  def scale_and_filter(self, experiments, reflections, params, scaling_params, filtering_params):
    """Write the behaviour of the program as functions and classes outside run()"""

    results = AnalysisResults()

    for counter in range(1, params.filtering.max_cycles+1):
      # First run scaling, followed by filtering
      scaling_script = self.scaling_script(scaling_params, experiments, reflections)
      scaling_script.run(save_data=False)
      # Log initial expids here, need to do after dataset selection in scaling
      # but before any filtering
      if counter == 1:
        results.initial_expids_and_image_ranges = [(exp.identifier, exp.scan.get_image_range()) \
          if exp.scan else None for exp in experiments]
      filter_script = self.filtering_script(filtering_params, scaling_script.experiments,
        scaling_script.reflections)
      filter_script.run()

      # Log results from scaling and filtering
      results = self.log_cycle_results(results, scaling_script, filter_script)
      logger.info("Cycle %s of filtering, n_reflections removed this cycle: %s",
        counter, results.get_last_cycle_results()['n_removed'])

      # Reset dataset inclusion/exclusion to avoid errors for repeated scaling.
      scaling_params.dataset_selection.use_datasets = None
      scaling_params.dataset_selection.exclude_datasets = None

      # Test termination conditions
      latest_results = results.get_last_cycle_results()
      if latest_results['n_removed'] == 0:
        logger.info("Finishing scale and filtering as no data removed in this cycle.")
        scaling_script.output()
        results.finish(termination_reason='no_more_removed')
        break

      if latest_results['cumul_percent_removed'] > params.filtering.max_percent_removed:
        logger.info("Finishing scale and filtering as have now removed more than the limit.")
        results = self._run_final_scale(
          scaling_params, experiments, reflections, results)
        results.finish(termination_reason='max_percent_removed')
        break

      if params.filtering.min_completeness:
        if latest_results['merging_stats']['completeness'] < params.filtering.min_completeness:
          logger.info("Finishing scale and filtering as completeness now below cutoff.")
          results = self._run_final_scale(
            scaling_params, experiments, reflections, results)
          results.finish(termination_reason='below_completeness_limit')
          break

      if counter == params.filtering.max_cycles:
        logger.info("Finishing as reached max number of cycles.")
        results = self._run_final_scale(
          scaling_params, experiments, reflections, results)
        results.finish(termination_reason='max_cycles')
        break

    # Print summary of results
    logger.info("\nSummary of data removed:")
    for i, res in enumerate(results.get_cycle_results()):
      logger.info("Cycle number: %s", i+1)
      if 'image_ranges_removed' in res:
        if res['image_ranges_removed']:
          logger.info("  Removed image ranges: \n    %s", '\n    '.join(
            str(t[0]) + ', dataset '+str(t[1]) for t in res['image_ranges_removed']))
      else:
        if res['removed_datasets']:
          logger.info("  Removed datasets: %s", res['removed_datasets'])
      logger.info("  cumulative %% of reflections removed: %s",
        res['cumul_percent_removed'])

    return results

  def _run_final_scale(self, scaling_params, experiments, reflections, results):
    scaling_script = self.scaling_script(scaling_params, experiments, reflections)
    scaling_script.run(save_data=True)
    results.add_final_stats(scaling_script.merging_statistics_result)
    return results

  @staticmethod
  def log_cycle_results(results, scaling_script, filter_script):
    """Log results from the scripts for this cycle and add to the results."""
    cycle_results = {'merging_stats' : scaling_script.merging_statistics_result}

    if not results.get_cycle_results():
      results.initial_n_reflections = scaling_script.scaled_miller_array.size()

    cycle_results['delta_cc_half_values'] = filter_script.results_summary[
      'per_dataset_delta_cc_half_values']['delta_cc_half_values']
    removal_summary = filter_script.results_summary['dataset_removal']
    if removal_summary['mode'] == 'image_group':
      cycle_results['image_ranges_removed'] = removal_summary['image_ranges_removed']
    cycle_results['removed_datasets'] = removal_summary['experiments_fully_removed']

    cycle_results['n_removed'] = filter_script.results_summary[
    'dataset_removal']['n_reflections_removed']

    n_removed = sum([res['n_removed'] for res in results.get_cycle_results()]) + \
      cycle_results['n_removed']
    percent_removed = n_removed/results.initial_n_reflections * 100
    cycle_results['cumul_percent_removed'] = percent_removed

    results.add_cycle_results(cycle_results)
    return results

def make_plots(analysis_results, experiments, params):
  """Make the three plots."""
  make_histogram_plots(analysis_results, params)
  make_merging_stats_plots(analysis_results, params)
  if params.delta_cc_half.mode == 'image_group':
    make_reduction_plots(analysis_results, experiments, params)

def make_histogram_plots(analysis_results, params):
  """Make and save the histogram plots."""
  cycle_results = analysis_results.get_cycle_results()
  delta_cc_half_lists = [res['delta_cc_half_values'] for res in cycle_results]

  if not delta_cc_half_lists:
    return
  n = len(delta_cc_half_lists)

  n_rows = int(floor(n/2.0)) + (n % 2)
  grid = gridspec.GridSpec(n_rows, 2)
  plt.figure(figsize=(10, 5*n_rows))
  axs = [plt.subplot(grid[int(c//2), int(c % 2)]) for c in range(n)]

  color_list = ['#F44336', '#FFC107', '#FFEB3B', '#8BC34A', '#03A9F4', '#3F51B5', '#607D8B']
  colors = [color_list[i%7] for i in range(n)]
  ordinal = lambda n: "%d%s" % (n, "tsnrhtdd"[(n/10%10 != 1)*(n%10 < 4)*n%10::4])
  legends = [ordinal(i)+r' $\Delta$CC$_{1/2}$'+"\nanalysis" for i in range(1, n+1)]
  if 'image_ranges_removed' in cycle_results[0]:
    n_rej = [len(res['image_ranges_removed']) for res in cycle_results]
  else:
    n_rej = [len(res['removed_datasets']) for res in cycle_results]

  for i, (ax, deltas) in enumerate(
    zip(axs, delta_cc_half_lists)):
    counts, bins, _ = ax.hist(np.array(deltas), bins=40, color=colors[i], label=legends[i])
    bin_width = bins[1] - bins[0]
    ymax = max(4, min(max(counts), 25))
    ax.set_ylim([0, ymax])
    ax.set_xlabel(r'$\Delta$CC$_{1/2}$')
    ax.set_ylabel('No. of image groups')
    ax.set_alpha(0.2)
    # Now plot some arrows
    xmin, xmax, ymin, ymax = ax.axis()
    n = 0
    for count, b in zip(counts, bins):
      if n >= n_rej[i]:
        break
      if count > 0:
        headl = (ymax-ymin)/16.0
        ax.arrow(b+(bin_width)/2.0, ymax/2.0, 0, -ymax/4.0 + 1.1*headl,
          width=(xmax-xmin)/(5*80.0), head_width=(xmax-xmin)/20.0,
          head_length=headl, color='k')
        n += count
    ax.legend(shadow=True, facecolor='w', framealpha=1, loc=2)
    #ax.background(False)
  plt.tight_layout()
  logger.info("Saving histogram plots to %s", params.output.plots.histogram)
  plt.savefig(params.output.plots.histogram)

def make_merging_stats_plots(analysis_results, params):
  """Make and save the merging stats plots."""
  merging_stats = analysis_results.get_merging_stats()
  plt.figure(figsize=(14, 8))
  grid = gridspec.GridSpec(2, 3)
  axs = [plt.subplot(grid[int(c//3), int(c % 3)]) for c in range(3)]
  color_list = ['#F44336', '#FFC107', '#FFEB3B', '#8BC34A', '#03A9F4', '#3F51B5', '#607D8B']
  colors = [color_list[i%7] for i in range(len(merging_stats))]
  colors[-1] = 'k'
  markers = ['s', 'o', 'v', '*']
  ordinal = lambda n: "%d%s" % (n, "tsnrhtdd"[(n/10%10 != 1)*(n%10 < 4)*n%10::4])
  legends = ['initial scale']
  if len(merging_stats) > 2:
    legends += [ordinal(i)+" rescale" for i in range(1, len(merging_stats)-1)] + ['final rescale']
  elif len(merging_stats) == 2:
    legends += ['final rescale']
  for c, stats in enumerate(merging_stats):
    x = [uctbx.d_as_d_star_sq(d) for d in stats['d_min']]
    axs[0].scatter(x, stats['ccs'], color=colors[c],
      marker=markers[c%4], label=legends[c])
    axs[0].plot(x, stats['ccs'], color=colors[c], label=None)
    axs[1].scatter(x, stats['rmerge'], color=colors[c],
      marker=markers[c%4], label=legends[c])
    axs[1].plot(x, stats['rmerge'], color=colors[c], label=None)
    axs[2].scatter(x, stats['rpim'], color=colors[c],
      marker=markers[c%4], label=legends[c])
    axs[2].plot(x, stats['rpim'], color=colors[c], label=None)
  xticks = axs[0].get_xticks()
  new_labels = [str(round(1/x**0.5, 2)) if x != 0 else str(0) for x in xticks]
  ylabels = ['CC$_{1/2}$', r'R$_{merge}$', r'R$_{pim}$']
  ylims = [1.05, min(1.5, max(axs[1].get_ylim())), min(1.5, max(axs[2].get_ylim()))]
  legend_locs = [3, 2, 2]
  for i, ax in enumerate(axs):
    ax.set_xticklabels(new_labels)
    ax.set_xlabel(r'd ($\AA$)')
    ax.grid()
    ax.set_ylabel(ylabels[i])
    ax.set_ylim([0.0, ylims[i]])
    ax.legend(shadow=True, facecolor='w', framealpha=1, loc=legend_locs[i])
  ax3 = plt.subplot(grid[1, :])
  cell_text = []
  columns = ['CC$_{1/2}$', r'R$_{merge}$', r'R$_{pim}$', '<I/sigma>',
    'completeness', 'n observations']
  for c, stats in enumerate(merging_stats):
    cell_text.append(['%.3f' % v if k != 'n_obs' else str(int(v)) \
      for k, v in stats['overall'].iteritems()])
  ax3.axis('off')
  ax3.table(cellText=cell_text, rowLabels=legends, colLabels=columns, loc='center')
  logger.info("Saving merging statistics plots to %s", params.output.plots.merging_stats)
  plt.savefig(params.output.plots.merging_stats)

def make_reduction_plots(analysis_results, experiments, params):
  """Make a chart showing excluded image ranges."""
  valid_image_ranges = get_valid_image_ranges(experiments)
  expids = [exp.identifier for exp in experiments]
  # first plot all initial image ranges
  x1 = range(len(analysis_results.initial_expids_and_image_ranges))
  y1_tops = []
  y1_bottoms = []
  initial_expids = []
  for expid_and_img in analysis_results.initial_expids_and_image_ranges:
    y1_tops.append(expid_and_img[1][1] - expid_and_img[1][0] + 1)
    y1_bottoms.append(expid_and_img[1][0] - 1)
    initial_expids.append(expid_and_img[0])
  x2 = []
  y2_tops = []
  y2_bottoms = []
  for (valid, expid) in zip(valid_image_ranges, expids):
    loc = initial_expids.index(expid)
    for t in valid:
      y2_tops.append(t[1]-t[0]+1)
      y2_bottoms.append(t[0]-1)
      x2.append(loc)
  plt.figure(figsize=(7, 7))
  plt.bar(x1, y1_tops, bottom=y1_bottoms, width=1, color='k', alpha=0.5, label="Before filtering")
  plt.bar(x2, y2_tops, bottom=y2_bottoms, width=1, color='b', alpha=0.5, label="After filtering")
  plt.xlabel("Dataset")
  plt.ylabel("Image number")
  plt.legend(loc=1)
  logger.info("Saving image range plots to %s", params.output.plots.image_ranges)
  plt.savefig(params.output.plots.image_ranges)


def run(args=None, phil=phil_scope):
  """Run the scale and filter script."""
  import dials.util.log
  from dials.util.options import OptionParser
  from dials.util.options import flatten_reflections
  from dials.util.options import flatten_experiments

  usage = "dials.scale_and_filter [options] integrated_experiments.json integrated.pickle"

  parser = OptionParser(usage=usage, phil=phil, read_reflections=True,
    read_experiments=True, check_format=False, epilog=__doc__)
  params, _ = parser.parse_args(args=args, show_diff_phil=False)

  scaling_params = params.scale
  filtering_params = params.delta_cc_half

  dials.util.log.config(info=params.output.log, debug=params.output.debug.log)

  diff_phil = parser.diff_phil.as_str()
  if diff_phil:
    logger.info('The following parameters have been modified:\n%s', diff_phil)

  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)

  scale_and_filter = ScaleAndFilter(ScalingScript, FilterScript)
  analysis_results = scale_and_filter.scale_and_filter(
    experiments, reflections, params, scaling_params, filtering_params)

  with open(params.output.analysis_results, 'w') as f:
    json.dump(analysis_results.to_dict(), f, indent=2)
  make_plots(analysis_results, experiments, params)

if __name__ == '__main__':
  with dials.util.show_mail_on_error():
    run()
