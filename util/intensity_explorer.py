#!/usr/bin/env python
#-*- coding:utf-8 -*-

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

Example:
  $ dials.python intensity_explorer.py <unmerged MTZ file>
"""

# TODO Docstrings are in Google-ish format — move to Sphinx-ish.
# TODO Docstrings need updating since new reflection table functionality.
# TODO Once ∃ a dials tool for (unmerged MTZ) –> (exp list, refl table), use it

import os
import sys
import math
from scipy import stats as ss
from cctbx import miller
from iotbx import mtz
from dials.array_family import flex
from matplotlib import colors, pyplot as plt


class IntensityDist(object):
  # FIXME Update docstring.
  """
  Store intensity data and generate normal order statistic medians.
  
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
    Imeans (cctbx_array_family_flex_ext.double):  Weighted means of symmetry-
                                                    equivalent reflection
                                                    intensities.
    sigImeans
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
    ind_unique (set):     Set of observed symmetry-inequivalent Miller indices.
    kept_singles (bool):  Indicates whether multiplicity-1 reflections were
                            retained.
                            Defaults to False.
    outfile (str):        File root for generated plots.
                            Defaults to MTZ input file root.
  """

  def __init__(self, rtable=None, elist=None, filename=None, outfile=None,
               keep_singles=False, uncertainty='sigma'):
    # FIXME Update docstring
    """
    Generate z-scores and normal probability plot from an unmerged MTZ file
    
    Args:
      filename (str):       Unmerged MTZ input file.
      outfile (str):        File root for output PNG plots.  If None, the root
                              of the input filename is used.
                              Defaults to None.
      keep_singles (bool):  Choose whether to keep multiplicity-1 reflections.
                              Defaults to False.
      uncertainty (str):    Measure of spread to use in normalising the
                              z-scores, i.e. z = (I - <I>) / uncertainty.
        Possible values for uncertainty:
        'sigma':    Use measured sigma values;
        'stddev':   Use sample standard deviations calculated as square-root of
          unbiased weighted sample variances of symmetry-equivalent reflection
          intensities;
        'sigImean': Use standard deviation on the weighted mean intensities.
          Mathematically meaningless, this is just for debugging.
        Defaults to 'sigma'.
    """

    if outfile:
      self.outfile = os.path.splitext(os.path.basename(outfile))[0]
    elif filename:
      self.outfile = os.path.splitext(os.path.basename(filename))[0]
    if rtable and elist:
      # FIXME Handle ExperimentList with length > 1 properly.
      self.rtable = rtable
      self.sg_type = elist.crystals()[0].get_space_group().type()
    elif filename:
      self.rtable = self._data_from_unmerged_mtz(filename)
    else:
      raise AttributeError(
          'Please specify either rtable and elist, or an input filename.')
    self.ind_unique = self._determine_multiplicity()
    if not keep_singles:
      self.ind_unique = self._discard_singletons()
    self._mean_error_stddev()
    self._make_z(uncertainty)
    self._probplot_data()

  def _data_from_unmerged_mtz(self, filename):
    """
    Produce a minimal reflection table from an MTZ file.

    The returned reflection table will not contain all the standard
    columns, only those that are necessary for this script.  Sorry!

    :param filename: :type str: Name of an unmerged MTZ input file.
    :return: :type dials.array_family_flex_ext.reflection_table: A
      reflection table object, containing only the columns
      * ``miller_index``
      * ``intensity.sum.value``
      * ``intensity.sum.variance``
      * ``xyzobs.px.value``
    """

    m = mtz.object(filename)  #Parse MTZ, with lots of useful methods.

    ind = m.extract_miller_indices()  #A flex array of Miller indices.

    cols = m.columns()  #Generates columns (augmented flex arrays).
    col_dict = { c.label() : c for c in cols }  #A dict of all the columns.
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

  def _data_from_pickle(self):
    # TODO Make this do stuff.
    pass

  def _determine_multiplicity(self):
    self.rtable['miller_index.asu'] = self.rtable['miller_index'].deep_copy()
    # Map Miller indices to asymmetric unit.
    # TODO Handle anomalous flag sensibly.  Currently assumes not anomalous.
    miller.map_to_asu(self.sg_type, False, self.rtable['miller_index.asu'])

    # Find unique Miller indices.
    ind_unique = set(self.rtable['miller_index.asu'])

    # Record the multiplicities.
    self.rtable['multiplicity'] = flex.int(self.rtable.size(), 0)
    for hkl in ind_unique:
      sel = (self.rtable['miller_index.asu'] == hkl).iselection()
      self.rtable['multiplicity'].set_selected(sel, sel.size())

    return ind_unique

  def _discard_singletons(self):
    # Drop multiplicity-1 data.
    sel = (self.rtable['multiplicity'] != 1).iselection()
    self.rtable = self.rtable.select(sel)
    ind_unique = set(self.rtable['miller_index.asu'])

    return ind_unique

  def _mean_error_stddev(self):
    # Calculate the weighted mean intensities.
    Imeans = flex.double(self.rtable.size(), 0)
    # Calculate the standard deviations from unbiased weighted variances.
    variances = flex.double(self.rtable.size(), 0)
    # For weighted averaging.
    weights = 1 / self.rtable['intensity.sum.variance']
    sum_weights = flex.double(self.rtable.size(), 0)
    sum_square_weights = flex.double(self.rtable.size(), 0)
    weighted_sum_square_residuals = flex.double(self.rtable.size(), 0)

    for hkl in self.ind_unique:
      sel = (self.rtable['miller_index'] == hkl).iselection()
      weight = weights.select(sel)
      I = self.rtable['intensity.sum.value'].select(sel)
      sum_weight = flex.sum(weight)
      try:
        Imean = flex.sum(weight * I) / sum_weight
      except:
        print(hkl, sel.size(), weight, I, sum_weight)
      weighted_sum_square_residual = flex.sum(weight * flex.pow2(I - Imean))

      sum_weights.set_selected(sel, sum_weight)
      Imeans.set_selected(sel, Imean)
      weighted_sum_square_residuals.set_selected(sel,
                                                 weighted_sum_square_residual)

    # Calculate the standard errors on the means.
    sigImeans = flex.sqrt(1 / sum_weights)
    # Treat variances differently for multiplicity-1 and higher multiplicity.
    sel = self.rtable.select(sel)['multiplicity'] != 1
    # For multiplicity-1 cases, just use the measured variance.
    variances.set_selected(~sel,
                           self.rtable['intensity.sum.variance'].select(~sel))
    # For higher multiplicity cases, calculate the unbiased weighted variance
    variances.set_selected(sel,
                           (weighted_sum_square_residuals.select(sel)
                            / (sum_weights
                               - sum_square_weights / sum_weights).select(sel))
                           )

    self.rtable['intensity.mean.value'] = Imeans
    self.rtable['intensity.mean.std_error'] = sigImeans
    self.rtable['intensity.mean.variance'] = variances

  def _make_z(self, uncertainty='sigma'):
    if uncertainty in ('stddev', 'sigImean'):
      uncertainty = {'stddev':
                       flex.sqrt(self.rtable['intensity.mean.variance']),
                     'sigImean':
                       self.rtable['intensity.mean.std_error']}[uncertainty]
    else:
      uncertainty = flex.sqrt(self.rtable['intensity.sum.variance'])

      z = (self.rtable['intensity.sum.value']
           - self.rtable['intensity.mean.value']) / uncertainty
      self.rtable['intensity.z_score'] = z

  def _probplot_data(self):
    order = flex.sort_permutation(self.rtable['intensity.z_score'])
    osm = flex.double(self.z.size(), 0)
    osm.set_selected(order, flex.double(ss.probplot(self.z, fit=False)[0]))
    self.rtable['intensity.order_statistic_medians'] = osm

  def plot_z_histogram(self):
    """Plot a histogram of the z-scores of the weighted mean intensities."""
    fig, ax = plt.subplots()

    ax.set_title(r'$z$ histogram')
    ax.set_xlabel(r'$z$')
    ax.set_ylabel(r'$N$')
    ax.hist(self.z, label='$z$', bins=100, range=(-5, 5))
    fig.savefig(self.outfile + '_zhistogram', transparent = True)
    plt.close()

  def _plot_symmetry_equivalents(self,
                                 overlay_mean=False, uncertainty=False
                                 ):
    """Really just a test plot.  Slow.  You probably don't want to use."""
    fig, ax = plt.subplots()

    for hkl in self.ind_unique:
      sel = (self.ind == hkl).iselection()

      if overlay_uncertainty == 'sigImean':
        yerr = self.sigImeans.select(sel)
      elif overlay_uncertainty:
        yerr = flex.sqrt(self.variances.select(sel))
      else:
        yerr = None

      plt.errorbar(
          sel,
          self.rtable['intensity.sum.value'].select(sel),
          yerr=flex.sqrt(self.rtable['intensity.sum.variance'].select(sel)),
          ls="--"
      )
      if overlay_mean:
        plt.errorbar(
            sel,
            self.Imeans.select(sel),
            yerr = yerr,
            ls = "-",
            color = "k",
            lw = .5
        )

    #fig.savefig(self.outfile + 'testfig')
    plt.show()

  def probplot(self, **kwargs):
    """Create a normal probability plot from the z-scores."""
    fig, ax = plt.subplots()

    ax.set_title('Normal probability plot')
    ax.set_xlabel('Order statistic medians, $m$')
    ax.set_ylabel(r'Ordered responses, $z$')
    ax.set_ylim(-10,10)
    ax.plot(self.osm, self.z, '.b', **kwargs)
    ax.plot([-5,5], [-5,5], '-g')

    fig.savefig(self.outfile + '_probplot', transparent = True)
    plt.close()

  def plot_z_vs_multiplicity(self, **kwargs):
    "Plot intensity z-scores versus multiplicity."
    fig, ax = plt.subplots()

    ax.set_title(r'$z$-scores versus multiplicity')
    ax.set_xlabel('Multiplicity')
    ax.set_ylabel(r'$z$')
    ax.plot(self.multis, self.z, '.', **kwargs)

    fig.savefig(self.outfile + '_z_vs_multiplicity', transparent = True)
    plt.close()

  def plot_z_map(self, minimum=0):
    """Plot a z-score heatmap of the detector.
    
    Beware, this is only meaningful if the data have a single geometry model.
    """
    sel = (flex.abs(self.z) >= minimum).iselection()
    x, y, z =  self.rtable.select(sel)['xyzobs.px.value'].parts()

    extreme = math.ceil(flex.max(flex.abs(self.z)))
    norm = colors.SymLogNorm(
        vmin = -extreme, vmax = extreme,
        linthresh = .02, linscale=1,
    )
    cmap_kws = {'cmap' : 'coolwarm_r', 'norm' : norm}

    fig, ax = plt.subplots()

    ax.set_title(r'$z$-scores versus detector position')
    ax.set_xlabel('Detector x position (pixels)')
    ax.set_ylabel('Detector y position (pixels)')
    ax.set_aspect('equal', 'box')
    ax.set_xlim(flex.min(x)-5, flex.max(x)+5)
    ax.set_ylim(flex.min(y)-5, flex.max(y)+5)
    det_map = ax.scatter(
        x,
        y,
        c = z,
        marker = ',',
        s = 0.5,
        **cmap_kws
    )
    cbar = fig.colorbar(det_map, ax=ax, **cmap_kws)
    cbar.set_label(r'$z$')

    fig.savefig(self.outfile + '_z_detector_map', transparent = True)
    plt.close()

  def plot_time_series(self, **kwargs):
    """Plot a crude time series of z-scores.
    
    Batch (frame) number is used as a proxy for time."""
    fig, ax = plt.subplots()

    ax.set_title(r'Time series of $z$-scores')
    ax.set_xlabel('Approximate chronology (frame number)')
    ax.set_ylabel(r'$z$')

    ax.plot(self.rtable['xyzobs.px.value'].parts()[2], self.z, '.', **kwargs)

    fig.savefig(self.outfile + '_z_time_series', transparent = True)
    plt.close()

  def plot_z_vs_IsigI(self, **kwargs):
    """Plot z-scores versus I/sigma."""
    fig, ax = plt.subplots()

    ax.set_title(r'$z$-scores versus spot intensity')
    ax.set_xlabel(r'$\bar{I}_\mathbf{h} / \sigma_\mathbf{h}$')
    ax.set_ylabel(r'$z$')
    ax.set_ylim(-10,10)
    ax.set_xscale('log')
    ax.plot(flex.abs(self.Imeans/self.sigImeans), self.z, '.', **kwargs)

    fig.savefig(self.outfile + '_z_vs_I_over_sigma', transparent = True)
    plt.close()

  def plot_z_vs_I(self, **kwargs):
    """Plot z-scores versus absolute intensity."""
    fig, ax = plt.subplots()

    ax.set_title(r'$z$-scores versus spot intensity')
    ax.set_xlabel(r'$\bar{I}_\mathbf{h}$')
    ax.set_ylabel(r'$z$')
    ax.set_ylim(-10,10)
    ax.set_xscale('log')
    ax.plot(flex.abs(self.Imeans), self.z, '.', **kwargs)

    fig.savefig(self.outfile + '_z_vs_I', transparent = True)
    plt.close()

if __name__ == "__main__":
  # TODO Handle multiple input MTZ files.
  # TODO Allow determination of output filename root.
  # Give an unmerged MTZ file as an argument:
  data = IntensityDist(filename=sys.argv[1])

  data.plot_z_histogram()
  data.probplot()
  data.plot_time_series()
  data.plot_z_map()
  data.plot_z_vs_multiplicity()
  data.plot_z_vs_I()
  data.plot_z_vs_IsigI()