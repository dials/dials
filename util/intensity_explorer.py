#!/usr/bin/env python
# -*- coding:utf-8 -*-

from __future__ import absolute_import, division, print_function

"""
Examine the distribution of diffraction spot intensities.

This module defines a class IntensityDist, with several methods for exploring
the distribution of measured spot intensities in an X-ray diffraction
experiment.  The user may wish to use this information to inform decisions
regarding the error model employed in analysing the data.  Data are passed in
as an unmerged MTZ file (see http://www.ccp4.ac.uk/html/mtzformat.html) and the
resulting IntensityDist instance contains the pertinent columns of data, along
with normal order statistic medians of the z-scores of the intensities, for
constructing a normal probability plot (See
https://www.itl.nist.gov/div898/handbook/eda/section3/normprpl.htm).

If called as a script, read data from an unmerged MTZ file; generate a
histogram and a normal probability plot of the z-scores of the intensity data,
along with plots of z as a function of batch number, of multiplicity, of
detector position, of measured multiplicity, of absolute intensity and of
I/sigma.
"""

# TODO Docstrings are in Google-ish format — move to Sphinx-ish.
# TODO Once ∃ a dials tool for (unmerged MTZ) –> (exp list, refl table), use it

from dials.array_family import flex


class IntensityDist(object):
  # TODO Finish updating docstring.
  """
  Store intensity data and generate associated normal probability plot data.

  :ivar outfile: Filename root for generated plots.  Default value: None
  :type outfile: str
  Attributes:
    ind
      (cctbx_array_family_flex_ext.miller_index): Miller indices.
    I (cctbx_array_family_flex_ext.double):       Measured intensity data.
    sigI (cctbx_array_family_flex_ext.double):    Measured intensity standard
                                                    deviations.
    x (cctbx_array_family_flex_ext.double):       Detector position, x (fast)
                                                    axis component.
    y (cctbx_array_family_flex_ext.double):       Detector position, y (slow)
                                                    axis component.
    frame (cctbx_array_family_flex_ext.int):      Batch (frame) number.
    multis (cctbx_array_family_flex_ext.int):     Measured multiplicity of
                                                    symmetry-equivalent spots.
    i_means (cctbx_array_family_flex_ext.double):  Weighted means of symmetry-
                                                    equivalent reflection
                                                    intensities.
    sig_i_means
      (cctbx_array_family_flex_ext.double):   Standard deviation on the
                                                weighted mean intensities.
    variances
      (cctbx_array_family_flex_ext.double):   Unbiased weighted sample
                                                variances of symmetry-
                                                equivalent reflection
                                                intensities.
    z (cctbx_array_family_flex_ext.double):   z-scores of weighted mean
                                                intensities.
    order
      (scitbx_array_family_flex_ext.size_t):  Index with which to sort the
                                                z-scores in ascending order.
                                                Useful for making a normal
                                                probability plot.
    osm (cctbx_array_family_flex_ext.double): Normal order statistic medians of
                                                the z-scores.
  """

  def __init__(self, rtable, elist, calculate_variances=False,
               keep_singles=False, uncertainty='sigma', outfile=None):
    # TODO Finish updating docstring
    """
    Generate z-scores and a normal probability plot from a DIALS
    reflection_table and a dxtbx ExperimentList, containing the observations
    and the corresponding experiments, respectively.

    :param rtable: A reflection table object, containing at least the columns
      * ``miller_index``
      * ``intensity.sum.value``
      * ``intensity.sum.variance``
      * ``xyzobs.px.value``
    :type rtable: dials.array_family_flex_ext.reflection_table
    :param elist: A corresponding experiment list.
    :type elist: dxtbx_model_ext.ExperimentList
    :param keep_singles: Choose whether to keep multiplicity-1
                          reflections.
                          Defaults to False.
    :type keep_singles: bool
    :param uncertainty: Measure of spread to use in normalising the
                        z-scores, i.e. z = (I - <I>) / uncertainty.
    :type uncertainty: str

    Possible values for uncertainty:
    * 'sigma':    Use measured sigma values;
    * 'stddev':   Use sample standard deviations calculated as
                  square-root of unbiased weighted sample variances
                  of symmetry-equivalent reflection intensities;
    Defaults to 'sigma'.
    """

    from cctbx import miller

    if not (rtable and elist):
      raise TypeError(
          "Must be called with a reflection table and experiment list.")

    # FIXME Handle ExperimentList with length > 1 properly.
    self.rtable = rtable
    self.sg_type = elist.crystals()[0].get_space_group().type()

    # Map Miller indices to asymmetric unit.
    self.rtable['miller_index.asu'] = self.rtable['miller_index'].deep_copy()
    # TODO Handle anomalous flag sensibly.  Currently assumes not anomalous.
    miller.map_to_asu(self.sg_type, False, self.rtable['miller_index.asu'])

    # Calculate normal probability plot data.
    self._multiplicity_mean_error_stddev(
      calculate_variances=calculate_variances, keep_singles=keep_singles)
    self._make_z(uncertainty)
    self._probplot_data()

    self.outfile = outfile

  def _multiplicity_mean_error_stddev(self, calculate_variances=False,
                                      keep_singles=False):
    u"""
    Calculate aggregate properties of grouped symmetry-equivalent reflections.

    Populate the reflection table of observations with the following
    properties:
      * ``multiplicity`` — Multiplicity of observations of a given reflection
      in the asymmetric unit;
      :type: `dials.array_family_flex_ext.int` array
      * ``intensity.mean.value`` — Mean of symmetry-equivalent reflections,
      weighted by measurement error;
      :type: `dials.array_family_flex_ext.double` array
      * ``intensity.mean.std_error`` — Standard error on the weighted mean;
      :type: `dials.array_family_flex_ext.double` array
      * (optional) ``intensity.mean.variance`` — variance of
      symmetry-equivalent reflections, weighted by measurement error;
      :type: `dials.array_family_flex_ext.double` array

    :param calculate_variances: Elect whether to calculate the weighted
    variances.  Defaults to False, to spare an expensive computation.
    :type calculate_variances: bool
    :param keep_singles: Choose whether to keep single-multiplicity
    reflections.
    :type keep_singles: bool
    """

    # Sort the reflection table for speedier iteration.
    self.rtable.sort('miller_index.asu')
    # Record the positions of any multiplicity-1 reflections.
    if not keep_singles:
      singles = flex.size_t()
    # Record the multiplicities.
    multiplicity = flex.int()
    # For weighted averaging.
    weights = 1 / self.rtable['intensity.sum.variance']
    sum_weights = flex.double()
    if calculate_variances:
      sum_square_weights = flex.double()
    # Calculate the weighted mean intensities.
    i_means = flex.double()
    # Calculate the standard deviations from unbiased weighted variances.
    variances = flex.double()

    # Iterate over the reflections, grouping by equivalent Miller index,
    # to calculate multiplicities, weighted mean intensities, etc..
    # Some time can be saved by only calculating variances if necessary.
    # Initial values:
    prev_index = None
    count = 1
    # One big loop through the entire reflection table:
    for j in range(self.rtable.size()):
      index = self.rtable['miller_index.asu'][j]
      weight = weights[j]
      # Aggregate within a symmetry-equivalent group of reflections:
      if index == prev_index:
        count += 1
        i_sum += weight * self.rtable['intensity.sum.value'][j]
        sum_weight += weight
        if calculate_variances:
          sum_square_weight += weight * weight
      # Record the aggregated values for the group:
      elif prev_index:
        if count == 1 and not keep_singles:
          singles.append(j-1)
        multiplicity.extend(flex.int(count, count))
        i_means.extend(flex.double(count, i_sum / sum_weight))
        sum_weights.extend(flex.double(count, sum_weight))
        if calculate_variances:
          sum_square_weights.extend(flex.double(count, sum_square_weight))
        # And reinitialise:
        prev_index = index
        count = 1
        i_sum = weight * self.rtable['intensity.sum.value'][j]
        sum_weight = weight
        if calculate_variances:
          sum_square_weight = weight * weight
      # Handle the first row:
      else:
        prev_index = self.rtable['miller_index.asu'][j]
        i_sum = weight * self.rtable['intensity.sum.value'][j]
        sum_weight = weight
        if calculate_variances:
          sum_square_weight = weight * weight
    # Record the aggregated values for the last group:
    if count == 1 and not keep_singles:
      singles.append(self.rtable.size() - 1)
    multiplicity.extend(flex.int(count, count))
    i_means.extend(flex.double(count, i_sum / sum_weight))
    sum_weights.extend(flex.double(count, sum_weight))
    if calculate_variances:
      sum_square_weights.extend(flex.double(count, sum_square_weight))

    # Discard singletons:
    if not keep_singles:
      singles_del = flex.bool(self.rtable.size(), True)
      singles_del.set_selected(singles, False)
      multiplicity, weights, sum_weights, i_means = [ a.select(singles_del)
        for a in multiplicity, weights, sum_weights, i_means ]
      self.rtable.del_selected(singles)
      if calculate_variances:
        sum_square_weights = sum_square_weights.select(singles_del)

    # Record the multiplicities in the reflection table.
    self.rtable['multiplicity'] = multiplicity
    # Record the weighted mean intensities in the reflection table.
    self.rtable['intensity.mean.value'] = i_means
    # Record the standard errors on the means in the reflection table.
    self.rtable['intensity.mean.std_error'] = flex.sqrt(1 / sum_weights)

    if calculate_variances:
      #Initialise values:
      prev_index = None
      for j in range(self.rtable.size()):
        index = self.rtable['miller_index.asu'][j]
        weight = weights[j]
        residual = self.rtable['intensity.sum.value'][j] - i_means[j]
        # Aggregate within a symmetry-equivalent group of reflections:
        if index == prev_index:
          count += 1
          weighted_sum_square_residual += weight * residual * residual
        # Record the aggregated value for the group:
        elif prev_index:
          # The weighted variance is undefined for multiplicity=1,
          # use the measured variance instead in this case.
          if count == 1:
            variances.append(self.rtable['intensity.sum.variance'][j-1])
          else:
            sum_weight = sum_weights[j-1]
            var_weight = 1 / (sum_weight -
                              sum_square_weights[j-1] / sum_weight)
            variances.extend(
                flex.double(count, weighted_sum_square_residual * var_weight))
          # Reinitialise:
          prev_index = index
          count = 1
          weighted_sum_square_residual = weight * residual * residual
        # Handle the first row:
        else:
          prev_index = self.rtable['miller_index.asu'][j]
          count = 1
          weighted_sum_square_residual = weight * residual * residual
      # Record the aggregated values for the last group:
      # The weighted variance is undefined for multiplicity=1,
      # use the measured variance instead in this case.
      if count == 1:
        variances.append(self.rtable['intensity.sum.variance'][-1])
      else:
        sum_weight = sum_weights[-1]
        var_weight = 1 / (sum_weight - sum_square_weights[-1] / sum_weight)
        variances.extend(
            flex.double(count, weighted_sum_square_residual * var_weight))
      # Record the variances in the reflection table.
      self.rtable['intensity.mean.variance'] = variances

  def _make_z(self, uncertainty='sigma'):
    u"""
    Generate reflection z-scores.

    Calculate z-scores from reflection intensities, weighted mean
    intensities and a chosen measure of uncertainty in the intensity
    measurement.

    :param uncertainty: Chosen measure of uncertainty.  Options are
      * ``stddev`` — standard deviation, as calculated from the unbiased
      weighted variance aggregated amongst all symmetry-equivalent reflections;
      * ``sigma`` — measurement error for individual reflections.
    :type uncertainty: str
    """

    uncertainty_name = {
      'sigma': 'intensity.sum.variance',
      'stddev': 'intensity.mean.variance'
    }[uncertainty]
    try:
      uncertainty_value = flex.sqrt(self.rtable[uncertainty_name])
    except KeyError:
      from logging import warn
      uncertainty_value = flex.sqrt(self.rtable['intensity.sum.variance'])
      warn(u"""Weighted variances haven't been calculated,
    be sure to specify calculate_variances=True to use them.
    Defaulting to measured σ values as a measure of uncertainty instead.""")

    z = (self.rtable['intensity.sum.value']
         - self.rtable['intensity.mean.value']) / uncertainty_value
    self.rtable['intensity.z_score'] = z

  def _probplot_data(self):
    """Generate the data for a normal probability plot of z-scores."""

    from scipy import stats as ss

    order = flex.sort_permutation(self.rtable['intensity.z_score'])
    osm = flex.double(self.rtable.size(), 0)
    osm.set_selected(order,
                     flex.double(ss.probplot(self.rtable['intensity.z_score'],
                                             fit=False)[0]))
    self.rtable['intensity.order_statistic_medians'] = osm

  def plot_z_histogram(self):
    """Plot a histogram of the z-scores of the weighted mean intensities."""

    from matplotlib import pyplot as plt

    fig, ax = plt.subplots()

    ax.set_title(r'$z$ histogram')
    ax.set_xlabel(r'$z$')
    ax.set_ylabel(r'$N$')
    ax.hist(self.rtable['intensity.z_score'],
            label='$z$', bins=100, range=(-5, 5))
    fig.savefig(self.outfile + '_zhistogram', transparent=True)
    plt.close()

  def _plot_symmetry_equivalents(self,
                                 overlay_mean=False,
                                 overlay_uncertainty=False):
    """Really just a test plot.  Slow.  You probably don't want to use."""

    from matplotlib import pyplot as plt

    ind_unique = set(self.rtable['miller_index.asu'])

    fig, ax = plt.subplots()

    for hkl in ind_unique:
      sel = (self.rtable['miller_index.asu'] == hkl).iselection()
      rtable_selected = self.rtable.select(sel)

      if overlay_uncertainty == 'sigImean':
        yerr = rtable_selected['intensity.mean.std_error']
      elif overlay_uncertainty:
        yerr = flex.sqrt(rtable_selected['intensity.mean.variance'])
      else:
        yerr = None

      plt.errorbar(sel,
                   rtable_selected['intensity.sum.value'],
                   yerr=flex.sqrt(rtable_selected['intensity.sum.variance']),
                   ls="--")
      if overlay_mean:
        plt.errorbar(sel,
                     rtable_selected['intensity.mean.value'],
                     yerr=yerr,
                     ls="-",
                     color="k",
                     lw=.5)

    #fig.savefig(self.outfile + 'testfig')
    #plt.close()
    plt.show()

  def probplot(self, **kwargs):
    """Create a normal probability plot from the z-scores."""

    from matplotlib import pyplot as plt

    fig, ax = plt.subplots()

    ax.set_title('Normal probability plot')
    ax.set_xlabel('Order statistic medians, $m$')
    ax.set_ylabel(r'Ordered responses, $z$')
    ax.set_ylim(-10, 10)
    ax.plot(self.rtable['intensity.order_statistic_medians'],
            self.rtable['intensity.z_score'],
            '.b', **kwargs)
    ax.plot([-5, 5], [-5, 5], '-g')

    fig.savefig(self.outfile + '_probplot', transparent=True)
    plt.close()

  def plot_z_vs_multiplicity(self, **kwargs):
    """Plot intensity z-scores versus multiplicity."""

    from matplotlib import pyplot as plt

    fig, ax = plt.subplots()

    ax.set_title(r'$z$-scores versus multiplicity')
    ax.set_xlabel('Multiplicity')
    ax.set_ylabel(r'$z$')
    ax.plot(self.rtable['multiplicity'], self.rtable['intensity.z_score'],
            '.', **kwargs)

    fig.savefig(self.outfile + '_z_vs_multiplicity', transparent=True)
    plt.close()

  def plot_z_map(self, minimum=0):
    """
    Plot a z-score heatmap of the detector.

    Beware, this is only meaningful if the data have a single geometry model.
    """

    import math
    from matplotlib import colors, pyplot as plt

    sel = (flex.abs(self.rtable['intensity.z_score']) >= minimum).iselection()
    x, y, z = self.rtable.select(sel)['xyzobs.px.value'].parts()

    extreme = math.ceil(flex.max(flex.abs(self.rtable['intensity.z_score'])))
    norm = colors.SymLogNorm(vmin=-extreme, vmax=extreme,
                             linthresh=.02, linscale=1,)
    cmap_kws = {'cmap': 'coolwarm_r', 'norm': norm}

    fig, ax = plt.subplots()

    ax.set_title(r'$z$-scores versus detector position')
    ax.set_xlabel('Detector x position (pixels)')
    ax.set_ylabel('Detector y position (pixels)')
    ax.set_aspect('equal', 'box')
    ax.set_xlim(flex.min(x) - 5, flex.max(x) + 5)
    ax.set_ylim(flex.min(y) - 5, flex.max(y) + 5)
    det_map = ax.scatter(x, y, c=z, marker=',', s=0.5, **cmap_kws)
    cbar = fig.colorbar(det_map, ax=ax, **cmap_kws)
    cbar.set_label(r'$z$')

    fig.savefig(self.outfile + '_z_detector_map', transparent=True)
    plt.close()

  def plot_time_series(self, **kwargs):
    """
    Plot a crude time series of z-scores.

    Batch (frame) number is used as a proxy for time."""

    from matplotlib import pyplot as plt

    fig, ax = plt.subplots()

    ax.set_title(r'Time series of $z$-scores')
    ax.set_xlabel('Approximate chronology (frame number)')
    ax.set_ylabel(r'$z$')

    ax.plot(self.rtable['xyzobs.px.value'].parts()[2],
            self.rtable['intensity.z_score'],
            '.', **kwargs)

    fig.savefig(self.outfile + '_z_time_series', transparent=True)
    plt.close()

  def plot_z_vs_IsigI(self, **kwargs):
    """Plot z-scores versus I/sigma."""

    from matplotlib import pyplot as plt

    fig, ax = plt.subplots()

    ax.set_title(r'$z$-scores versus spot intensity')
    ax.set_xlabel(r'$\bar{I}_\mathbf{h} / \sigma_\mathbf{h}$')
    ax.set_ylabel(r'$z$')
    ax.set_ylim(-10, 10)
    ax.set_xscale('log')
    ax.plot(flex.abs(self.rtable['intensity.mean.value']
                     / self.rtable['intensity.mean.std_error']),
            self.rtable['intensity.z_score'],
            '.', **kwargs)

    fig.savefig(self.outfile + '_z_vs_I_over_sigma', transparent=True)
    plt.close()

  def plot_z_vs_I(self, **kwargs):
    """Plot z-scores versus absolute intensity."""

    from matplotlib import pyplot as plt

    fig, ax = plt.subplots()

    ax.set_title(r'$z$-scores versus spot intensity')
    ax.set_xlabel(r'$\bar{I}_\mathbf{h}$')
    ax.set_ylabel(r'$z$')
    ax.set_ylim(-10, 10)
    ax.set_xscale('log')
    ax.plot(flex.abs(self.rtable['intensity.mean.value']),
            self.rtable['intensity.z_score'],
            '.', **kwargs)

    fig.savefig(self.outfile + '_z_vs_I', transparent=True)
    plt.close()


def data_from_unmerged_mtz(filename):
  """
  Produce a minimal reflection table from an MTZ file.

  The returned reflection table will not contain all the standard
  columns, only those that are necessary for the IntensityDist class.

  :param filename: Name of an unmerged MTZ input file.
  :type filename: str
  :return: A reflection table object, containing only the columns
    * ``miller_index``
    * ``intensity.sum.value``
    * ``intensity.sum.variance``
    * ``xyzobs.px.value``
  :rtype: dials.array_family_flex_ext.reflection_table
  """

  from iotbx import mtz

  m = mtz.object(filename)  # Parse MTZ, with lots of useful methods.

  ind = m.extract_miller_indices()  # A flex array of Miller indices.

  cols = m.columns()  # Generates columns (augmented flex arrays).
  col_dict = {c.label(): c for c in cols}  # A dict of all the columns.
  I, sigI, x, y = (
    col_dict[label].extract_values().as_double()
    for label in ('I', 'SIGI', 'XDET', 'YDET')
  )
  frame = col_dict['BATCH'].extract_values().as_double().iround().as_double()
  # Honestly flex?!  Oh well, for now, we have to go round the houses.

  rtable = flex.reflection_table()
  rtable['miller_index'] = ind
  rtable['intensity.sum.value'] = I
  rtable['intensity.sum.variance'] = flex.pow2(sigI)
  rtable['xyzobs.px.value'] = flex.vec3_double(x, y, frame)

  return rtable
  # TODO Make this produce a minimal experiment list too.


def data_from_pickle_and_json():
  pass


if __name__ == "__main__":
  import sys
  from dials.util.phil import ReflectionTableConverters
  from dxtbx.model.experiment_list import ExperimentListFactory

  # TODO Handle multiple input MTZ files/pairs of pickle & json files.
  # TODO Allow determination of output filename root.
  # Give a pickle and a json file as arguments:
  rtable = ReflectionTableConverters().from_string(sys.argv[1]).data
  elist = ExperimentListFactory.from_json_file(sys.argv[2])

  data = IntensityDist(rtable, elist, outfile='Test')
  data.plot_z_histogram()
  data.probplot()
  data.plot_time_series()
  data.plot_z_map()
  data.plot_z_vs_multiplicity()
  data.plot_z_vs_I()
  data.plot_z_vs_IsigI()
