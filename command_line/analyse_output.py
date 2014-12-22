#!/usr/bin/env python
#
# analyse_output.py
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
import matplotlib

# Offline backend
matplotlib.use("Agg")

from matplotlib import pylab

def ensure_directory(path):
  ''' Make the directory if not already there. '''
  from os import makedirs
  import errno
  try:
    makedirs(path)
  except OSError as e:
    if e.errno != errno.EEXIST:
      raise

def ensure_required(rlist, required):
  ''' Check which keys aren't present. '''
  not_present = []
  for k in required:
    if k not in rlist:
      not_present.append(k)
  if len(not_present) != 0:
    print " Skipping: following required fields not present:"
    for k in not_present:
      print "  %s" % k
    return False
  return True


class CentroidAnalyser(object):
  ''' Analyse the reflection centroids. '''

  def __init__(self, directory):
    ''' Setup the directory. '''
    from os.path import join

    # Set the directory
    self.directory = join(directory, "centroid")
    ensure_directory(self.directory)

    # Set the required fields
    self.required = [
      "intensity.sum.value",
      "intensity.sum.variance",
      "xyzcal.px",
      "xyzobs.px.value"
    ]

  def __call__(self, rlist):
    ''' Analyse the reflection centroids. '''
    from dials.util.command_line import Command

    # Check we have the required fields
    print "Analysing reflection centroids"
    if not ensure_required(rlist, self.required):
      return

    # Remove I_sigma <= 0
    selection = rlist['intensity.sum.variance'] <= 0
    if selection.count(True) > 0:
      rlist.del_selected(selection)
      print ' Removing %d reflections with variance <= 0' % \
        selection.count(True)

    # Select only integrated reflections
    Command.start(" Selecting only integated reflections")
    mask = rlist.get_flags(rlist.flags.integrated)
    if mask.count(True) > 0:
      threshold = 10
      rlist = rlist.select(mask)
      Command.end(" Selected %d integrated reflections" % len(rlist))
    else:
      # Select only those reflections used in refinement
      threshold = 0
      mask = rlist.get_flags(rlist.flags.used_in_refinement)
      rlist = rlist.select(mask)
      Command.end(" Selected %d refined reflections" % len(rlist))

    # Look at differences in calculated/observed position
    print " Analysing centroid differences with I/Sigma > %s" %threshold
    self.centroid_diff_hist(rlist, threshold)
    print " Analysing centroid differences in x/y with I/Sigma > %s" %threshold
    self.centroid_diff_xy(rlist, threshold)
    self.centroid_xy_xz_zy_residuals(rlist, threshold)
    print " Analysing centroid differences in z with I/Sigma > %s" %threshold
    self.centroid_diff_z(rlist, threshold)
    print " Analysing centroid differences vs phi with I/Sigma > %s" %threshold
    self.centroid_mean_diff_vs_phi(rlist, threshold)

  def centroid_diff_hist(self, rlist, threshold):
    ''' Analyse the correlations. '''
    from dials.array_family import flex
    from os.path import join
    I = rlist['intensity.sum.value']
    I_sig = flex.sqrt(rlist['intensity.sum.variance'])
    I_over_S = I / I_sig
    mask = I_over_S > threshold
    rlist = rlist.select(mask)
    assert(len(rlist) > 0)
    xc, yc, zc = rlist['xyzcal.px'].parts()
    xo, yo, zo = rlist['xyzobs.px.value'].parts()
    xd = xo - xc
    yd = yo - yc
    zd = zo - zc
    diff = flex.sqrt(xd*xd + yd*yd + zd*zd)
    pylab.title("Difference between observed and calculated")
    cax = pylab.hist(diff, bins=20)
    pylab.xlabel("Difference in position")
    pylab.ylabel("# reflections")
    pylab.savefig(join(self.directory, "centroid_diff_hist.png"))
    pylab.clf()

  def centroid_diff_xy(self, rlist, threshold):
    ''' Look at the centroid difference in x, y '''
    from dials.array_family import flex
    from os.path import join
    I = rlist['intensity.sum.value']
    I_sig = flex.sqrt(rlist['intensity.sum.variance'])
    I_over_S = I / I_sig
    mask = I_over_S > threshold
    rlist = rlist.select(mask)
    assert(len(rlist) > 0)
    xc, yc, zc = rlist['xyzcal.px'].parts()
    xo, yo, zo = rlist['xyzobs.px.value'].parts()
    xd = xo - xc
    yd = yo - yc
    pylab.title("Difference between observed and calculated in X")
    cax = pylab.hexbin(xc, yc, C=xd, gridsize=100)
    pylab.xlabel("x")
    pylab.ylabel("y")
    cbar = pylab.colorbar(cax)
    cbar.ax.set_ylabel("Difference in x position")
    pylab.savefig(join(self.directory, "centroid_diff_x.png"))
    pylab.clf()
    pylab.title("Difference between observed and calculated in Y")
    cax = pylab.hexbin(xc, yc, C=yd, gridsize=100)
    pylab.xlabel("x")
    pylab.ylabel("y")
    cbar = pylab.colorbar(cax)
    cbar.ax.set_ylabel("Difference in y position")
    pylab.savefig(join(self.directory, "centroid_diff_y.png"))
    pylab.clf()

  def centroid_diff_z(self, rlist, threshold):
    ''' Look at the centroid difference in x, y '''
    from dials.array_family import flex
    from os.path import join
    I = rlist['intensity.sum.value']
    I_sig = flex.sqrt(rlist['intensity.sum.variance'])
    I_over_S = I / I_sig
    mask = I_over_S > threshold
    rlist = rlist.select(mask)
    assert(len(rlist) > 0)
    xc, yc, zc = rlist['xyzcal.px'].parts()
    xo, yo, zo = rlist['xyzobs.px.value'].parts()
    zd = zo - zc
    pylab.title("Difference between observed and calculated in Z")
    cax = pylab.hexbin(zc, zd, gridsize=100)
    pylab.xlabel("z")
    pylab.ylabel("Difference in z position")
    cbar = pylab.colorbar(cax)
    cbar.ax.set_ylabel("# Reflections")
    pylab.savefig(join(self.directory, "centroid_diff_z.png"))
    pylab.clf()

  def centroid_mean_diff_vs_phi(self, rlist, threshold):
    from dials.array_family import flex
    from os.path import join
    import math
    I = rlist['intensity.sum.value']
    I_sig = flex.sqrt(rlist['intensity.sum.variance'])
    I_over_S = I / I_sig
    mask = I_over_S > threshold
    rlist = rlist.select(mask)
    assert(len(rlist) > 0)

    xc, yc, zc = rlist['xyzcal.mm'].parts()
    xo, yo, zo = rlist['xyzobs.mm.value'].parts()

    dx = xc - xo
    dy = yc - yo
    dphi = zc - zo

    mean_residuals_x = []
    mean_residuals_y = []
    mean_residuals_phi = []
    frame = []
    phi_obs_deg = (180/math.pi) * zo
    phi = []

    for i_phi in range(int(math.floor(flex.min(phi_obs_deg))),
                   int(math.ceil(flex.max(phi_obs_deg)))):
      sel = (phi_obs_deg >= i_phi) & (phi_obs_deg < (i_phi+1))
      if sel.count(True) == 0:
        continue
      mean_residuals_x.append(flex.mean(dx.select(sel)))
      mean_residuals_y.append(flex.mean(dy.select(sel)))
      mean_residuals_phi.append(flex.mean(dphi.select(sel)))
      phi.append(i_phi)

    from matplotlib import pyplot
    fig = pyplot.figure()
    ax = fig.add_subplot(311)
    fig.subplots_adjust(hspace=0.5)
    pyplot.axhline(0, color='grey')
    ax.scatter(phi, mean_residuals_x)
    ax.set_xlabel('phi (deg)')
    ax.set_ylabel('mean $\Delta$ x (mm)')
    ax = fig.add_subplot(312)
    pyplot.axhline(0, color='grey')
    ax.scatter(phi, mean_residuals_y)
    ax.set_xlabel('phi (deg)')
    ax.set_ylabel('mean $\Delta$ y (mm)')
    ax = fig.add_subplot(313)
    pyplot.axhline(0, color='grey')
    ax.scatter(phi, mean_residuals_phi)
    ax.set_xlabel('phi (deg)')
    ax.set_ylabel('mean $\Delta$ phi (deg)')
    pyplot.savefig(join(self.directory, "centroid_mean_diff_vs_phi.png"))
    pyplot.clf()

  def centroid_xy_xz_zy_residuals(self, rlist, threshold):
    from dials.array_family import flex
    from os.path import join
    import math
    I = rlist['intensity.sum.value']
    I_sig = flex.sqrt(rlist['intensity.sum.variance'])
    I_over_S = I / I_sig
    mask = I_over_S > threshold
    rlist = rlist.select(mask)
    assert(len(rlist) > 0)

    xc, yc, zc = rlist['xyzcal.mm'].parts()
    xo, yo, zo = rlist['xyzobs.mm.value'].parts()

    dx = xc - xo
    dy = yc - yo
    dphi = zc - zo

    panel_ids = rlist['panel']
    crystal_ids = rlist['id']
    n_crystals = flex.max(crystal_ids) + 1
    n_panels = flex.max(panel_ids) + 1
    n_cols = int(math.floor(math.sqrt(n_panels)))
    n_rows = int(math.ceil(n_panels / n_cols))

    from matplotlib import pyplot

    for i_crystal in range(n_crystals):
      if n_crystals > 1:
        suffix = '_%i' %i_crystal
      else:
        suffix = ''
      crystal_sel = (crystal_ids == i_crystal)
      fig_xy, axes_xy = pyplot.subplots(
        n_rows, n_cols, sharex=True, sharey=True)
      fig_xz, axes_xz = pyplot.subplots(
        n_rows, n_cols, sharex=True, sharey=True)
      fig_zy, axes_zy = pyplot.subplots(
        n_rows, n_cols, sharex=True, sharey=True)

      if n_panels == 1:
        axes_xy = [[axes_xy]]
        axes_xz = [[axes_xz]]
        axes_zy = [[axes_zy]]

      i_panel = 0
      for i_row in range(n_rows):
        for i_col in range(n_cols):
          panel_sel = (panel_ids == i_panel)

          axes_xy[i_row][i_col].axhline(0, color='grey')
          axes_xy[i_row][i_col].axvline(0, color='grey')
          ax_xy = axes_xy[i_row][i_col].scatter(
            dx.select(panel_sel & crystal_sel).as_numpy_array(),
            dy.select(panel_sel & crystal_sel).as_numpy_array(),
            c='b', alpha=0.3, label='Panel %d' %i_panel)
          axes_xy[i_row][i_col].axes.set_aspect('equal')

          axes_xz[i_row][i_col].axhline(0, color='grey')
          axes_xz[i_row][i_col].axvline(0, color='grey')
          ax_xz = axes_xz[i_row][i_col].scatter(
            dx.select(panel_sel & crystal_sel).as_numpy_array(),
            dphi.select(panel_sel & crystal_sel).as_numpy_array(),
            c='b', alpha=0.3, label='Panel %d' %i_panel)

          axes_zy[i_row][i_col].axhline(0, color='grey')
          axes_zy[i_row][i_col].axvline(0, color='grey')
          ax_zy = axes_zy[i_row][i_col].scatter(
            dphi.select(panel_sel & crystal_sel).as_numpy_array(),
            dy.select(panel_sel & crystal_sel).as_numpy_array(),
            c='b', alpha=0.3, label='Panel %d' %i_panel)

          if n_panels > 1:
            axes_xy[i_row][i_col].set_title('Panel %d' %i_panel)
            axes_xz[i_row][i_col].set_title('Panel %d' %i_panel)
            axes_zy[i_row][i_col].set_title('Panel %d' %i_panel)

          if (i_row+1) == n_rows:
            axes_xy[i_row][i_col].set_xlabel('$\Delta$(x) (mm)')
            axes_xz[i_row][i_col].set_xlabel('$\Delta$(x) (mm)')
            axes_zy[i_row][i_col].set_xlabel('$\Delta$(phi) (deg)')
          if i_col == 0:
            axes_xy[i_row][i_col].set_ylabel('$\Delta$(y) (mm)')
            axes_xz[i_row][i_col].set_ylabel('$\Delta$(phi) (deg)')
            axes_zy[i_row][i_col].set_ylabel('$\Delta$(y) (mm)')
          i_panel += 1

      for fig in (fig_xy, fig_xz, fig_zy):
        fig.set_size_inches([n_cols * i for i in fig.get_size_inches()])
      fig_xy.savefig(
        join(self.directory, "centroid_xy_residuals%s.png" %suffix))
      fig_xz.savefig(
        join(self.directory, "centroid_xz_residuals%s.png" %suffix))
      fig_zy.savefig(
        join(self.directory, "centroid_zy_residuals%s.png" %suffix))
      pyplot.clf()

class BackgroundAnalyser(object):
  ''' Analyse the background. '''

  def __init__(self, directory):
    ''' Setup the directory. '''
    from os.path import join

    # Set the directory
    self.directory = join(directory, "background")
    ensure_directory(self.directory)

    # Set the required fields
    self.required = [
      "background.mse",
      "background.mean",
      "intensity.sum.value",
      "intensity.sum.variance",
      "xyzcal.px",
    ]

  def __call__(self, rlist):
    ''' Analyse the relfection background. '''
    from dials.util.command_line import Command

    # Check we have the required fields
    print "Analysing reflection backgrounds"
    if not ensure_required(rlist, self.required):
      return

    selection = rlist['intensity.sum.variance'] <= 0
    if selection.count(True) > 0:
      rlist.del_selected(selection)
      print ' Removing %d reflections with variance <= 0' % \
        selection.count(True)

    selection = rlist['background.mse'] < 0
    if selection.count(True) > 0:
      rlist.del_selected(selection)
      print ' Removing %d reflections with negative background model RMSD' % \
        selection.count(True)

    selection = rlist['background.mean'] <= 0
    if selection.count(True) > 0:
      rlist.del_selected(selection)
      print ' Removing %d reflections with mean background <= 0' % \
        selection.count(True)

    # Select only integrated reflections
    Command.start(" Selecting only integated reflections")
    mask = rlist.get_flags(rlist.flags.integrated)
    if mask.count(True) == 0:
      return

    rlist = rlist.select(mask)
    Command.end(" Selected %d integrated reflections" % len(rlist))

    # Look at distribution of I/Sigma
    print " Analysing distribution of background mean"
    self.mean_hist(rlist)
    print " Analysing distribution of background mean vs XY"
    self.mean_vs_xy(rlist)
    print " Analysing distribution of background mean vs z"
    self.mean_vs_z(rlist)
    print " Analysing distribution of background mean vs I/Sigma"
    self.mean_vs_ios(rlist)
    print " Analysing distribution of background CVRMSD"
    self.rmsd_hist(rlist)
    print " Analysing distribution of background CVRMSD vs XY"
    self.rmsd_vs_xy(rlist)
    print " Analysing distribution of background CVRMSD vs z"
    self.rmsd_vs_z(rlist)
    print " Analysing distribution of background CVRMSD vs I/Sigma"
    self.rmsd_vs_ios(rlist)

  def mean_hist(self, rlist):
    ''' Analyse the background RMSD. '''
    from os.path import join
    MEAN = rlist['background.mean']
    pylab.title("Background Model mean histogram")
    pylab.hist(MEAN, bins=20)
    pylab.xlabel("mean")
    pylab.ylabel("# reflections")
    pylab.savefig(join(self.directory, "background_model_mean_hist"))
    pylab.clf()

  def mean_vs_xy(self, rlist):
    ''' Plot I/Sigma vs X/Y '''
    from os.path import join
    MEAN = rlist['background.mean']
    x, y, z = rlist['xyzcal.px'].parts()
    pylab.title("Distribution of Background Model mean vs X/Y")
    cax = pylab.hexbin(x, y, C=MEAN, gridsize=100)
    pylab.xlabel("x")
    pylab.ylabel("y")
    cbar = pylab.colorbar(cax)
    cbar.ax.set_ylabel("Background Model mean")
    pylab.savefig(join(self.directory, "background_model_mean_vs_xy.png"))
    pylab.clf()

  def mean_vs_z(self, rlist):
    ''' Plot I/Sigma vs Z. '''
    from os.path import join
    MEAN = rlist['background.mean']
    x, y, z = rlist['xyzcal.px'].parts()
    pylab.title("Distribution of Background Model mean vs Z")
    cax = pylab.hexbin(z, MEAN, gridsize=100)
    pylab.xlabel("z")
    pylab.ylabel("Background Model mean")
    cbar = pylab.colorbar(cax)
    cbar.ax.set_ylabel("# reflections")
    pylab.savefig(join(self.directory, "background_model_mean_vs_z.png"))
    pylab.clf()

  def mean_vs_ios(self, rlist):
    ''' Analyse the correlations. '''
    from dials.array_family import flex
    from os.path import join
    MEAN = rlist['background.mean']
    I = rlist['intensity.sum.value']
    I_sig = flex.sqrt(rlist['intensity.sum.variance'])
    I_over_S = I / I_sig
    mask = I_over_S > 0.1
    I_over_S = I_over_S.select(mask)
    MEAN = MEAN.select(mask)
    pylab.title("Background Model mean vs Log I/Sigma")
    cax = pylab.hexbin(flex.log(I_over_S), MEAN, gridsize=100)
    cbar = pylab.colorbar(cax)
    pylab.xlabel("Log I/Sigma")
    pylab.ylabel("Background Model mean")
    cbar.ax.set_ylabel("# reflections")
    pylab.savefig(join(self.directory, "background_model_mean_vs_ios.png"))
    pylab.clf()

  def rmsd_hist(self, rlist):
    ''' Analyse the background RMSD. '''
    from dials.array_family import flex
    from os.path import join
    RMSD = flex.sqrt(rlist['background.mse'])
    MEAN = rlist['background.mean']
    RMSD = RMSD / MEAN
    pylab.title("Background Model mean histogram")
    pylab.hist(RMSD, bins=20)
    pylab.xlabel("mean")
    pylab.ylabel("# reflections")
    pylab.savefig(join(self.directory, "background_model_cvrmsd_hist"))
    pylab.clf()

  def rmsd_vs_xy(self, rlist):
    ''' Plot I/Sigma vs X/Y '''
    from dials.array_family import flex
    from os.path import join
    RMSD = flex.sqrt(rlist['background.mse'])
    MEAN = rlist['background.mean']
    RMSD = RMSD / MEAN
    x, y, z = rlist['xyzcal.px'].parts()
    pylab.title("Distribution of Background Model CVRMSD vs X/Y")
    cax = pylab.hexbin(x, y, C=RMSD, gridsize=100)
    pylab.xlabel("x")
    pylab.ylabel("y")
    cbar = pylab.colorbar(cax)
    cbar.ax.set_ylabel("Background Model CVRMSD")
    pylab.savefig(join(self.directory, "background_model_cvrmsd_vs_xy.png"))
    pylab.clf()

  def rmsd_vs_z(self, rlist):
    ''' Plot I/Sigma vs Z. '''
    from dials.array_family import flex
    from os.path import join
    RMSD = flex.sqrt(rlist['background.mse'])
    MEAN = rlist['background.mean']
    RMSD = RMSD / MEAN
    x, y, z = rlist['xyzcal.px'].parts()
    pylab.title("Distribution of Background Model CVRMSD vs Z")
    cax = pylab.hexbin(z, RMSD, gridsize=100)
    pylab.xlabel("z")
    pylab.ylabel("Background Model CVRMSD")
    cbar = pylab.colorbar(cax)
    cbar.ax.set_ylabel("# reflections")
    pylab.savefig(join(self.directory, "background_model_cvrmsd_vs_z.png"))
    pylab.clf()

  def rmsd_vs_ios(self, rlist):
    ''' Analyse the correlations. '''
    from dials.array_family import flex
    from os.path import join
    RMSD = flex.sqrt(rlist['background.mse'])
    MEAN = rlist['background.mean']
    RMSD = RMSD / MEAN
    I = rlist['intensity.sum.value']
    I_sig = flex.sqrt(rlist['intensity.sum.variance'])
    I_over_S = I / I_sig
    mask = I_over_S > 0.1
    I_over_S = I_over_S.select(mask)
    RMSD = RMSD.select(mask)
    pylab.title("Background Model CVRMSD vs Log I/Sigma")
    cax = pylab.hexbin(flex.log(I_over_S), RMSD, gridsize=100)
    cbar = pylab.colorbar(cax)
    pylab.xlabel("Log I/Sigma")
    pylab.ylabel("Background Model CVRMSD")
    cbar.ax.set_ylabel("# reflections")
    pylab.savefig(join(self.directory, "background_model_cvrmsd_vs_ios.png"))
    pylab.clf()


class IntensityAnalyser(object):
  ''' Analyse the intensities. '''

  def __init__(self, directory):
    ''' Set up the directory. '''
    from os.path import join

    # Set the directory
    self.directory = join(directory, "intensity")
    ensure_directory(self.directory)

    # Set the required fields
    self.required = [
      "intensity.sum.value",
      "intensity.sum.variance",
      "xyzcal.px",
    ]

  def __call__(self, rlist):
    ''' Analyse the reflection centroids. '''
    from dials.util.command_line import Command

    # FIXME Do the same and a comparison for intensity.prf

    # Check we have the required fields
    print "Analysing reflection intensities"
    if not ensure_required(rlist, self.required):
      return

    selection = rlist['intensity.sum.variance'] <= 0
    if selection.count(True) > 0:
      rlist.del_selected(selection)
      print ' Removing %d reflections with variance <= 0' % \
        selection.count(True)

    selection = rlist['intensity.sum.value'] <= 0
    if selection.count(True) > 0:
      rlist.del_selected(selection)
      print ' Removing %d reflections with intensity <= 0' % \
        selection.count(True)

    # Select only integrated reflections
    Command.start(" Selecting only integated reflections")
    mask = rlist.get_flags(rlist.flags.integrated)
    if mask.count(True) == 0:
      return

    rlist = rlist.select(mask)
    Command.end(" Selected %d integrated reflections" % len(rlist))

    # Look at distribution of I/Sigma
    print " Analysing distribution of I/Sigma"
    self.i_over_s_hist(rlist)
    print " Analysing distribution of I/Sigma vs xy"
    self.i_over_s_vs_xy(rlist)
    print " Analysing distribution of I/Sigma vs z"
    self.i_over_s_vs_z(rlist)
    print " Analysing number of background pixels used"
    self.num_background_hist(rlist)
    print " Analysing number of foreground pixels used"
    self.num_foreground_hist(rlist)

  def i_over_s_hist(self, rlist):
    ''' Analyse the correlations. '''
    from dials.array_family import flex
    from os.path import join
    I = rlist['intensity.sum.value']
    I_sig = flex.sqrt(rlist['intensity.sum.variance'])
    I_over_S = I / I_sig
    pylab.title("Log I/Sigma histogram")
    pylab.hist(flex.log(I_over_S), bins=20)
    pylab.xlabel("Log I/Sigma")
    pylab.ylabel("# reflections")
    pylab.savefig(join(self.directory, "ioversigma_hist"))
    pylab.clf()

  def i_over_s_vs_xy(self, rlist):
    ''' Plot I/Sigma vs X/Y '''
    from dials.array_family import flex
    from os.path import join
    I = rlist['intensity.sum.value']
    I_sig = flex.sqrt(rlist['intensity.sum.variance'])
    I_over_S = I / I_sig
    x, y, z = rlist['xyzcal.px'].parts()
    pylab.title("Distribution of I/Sigma vs X/Y")
    cax = pylab.hexbin(x, y, C=flex.log(I_over_S), gridsize=100)
    pylab.xlabel("x")
    pylab.ylabel("y")
    cbar = pylab.colorbar(cax)
    cbar.ax.set_ylabel("Log I/Sigma")
    pylab.savefig(join(self.directory, "ioversigma_vs_xy.png"))
    pylab.clf()

  def i_over_s_vs_z(self, rlist):
    ''' Plot I/Sigma vs Z. '''
    from dials.array_family import flex
    from os.path import join
    I = rlist['intensity.sum.value']
    I_sig = flex.sqrt(rlist['intensity.sum.variance'])
    I_over_S = I / I_sig
    x, y, z = rlist['xyzcal.px'].parts()
    pylab.title("Distribution of I/Sigma vs Z")
    cax = pylab.hexbin(z, flex.log(I_over_S), gridsize=100)
    pylab.xlabel("z")
    pylab.ylabel("Log I/Sigma")
    cbar = pylab.colorbar(cax)
    cbar.ax.set_ylabel("# reflections")
    pylab.savefig(join(self.directory, "ioversigma_vs_z.png"))
    pylab.clf()

  def num_background_hist(self, rlist):
    ''' Analyse the number of background pixels. '''
    from os.path import join
    if 'n_background' in rlist:
      N = rlist['n_background']
      pylab.title("Num Background Pixel Histogram")
      pylab.hist(N, bins=20)
      pylab.xlabel("Number of pixels")
      pylab.ylabel("# reflections")
      pylab.savefig(join(self.directory, "n_background_hist.png"))
      pylab.clf()

  def num_foreground_hist(self, rlist):
    ''' Analyse the number of foreground pixels. '''
    from os.path import join
    if 'n_foreground' in rlist:
      N = rlist['n_foreground']
      pylab.title("Num Foreground Pixel Histogram")
      pylab.hist(N, bins=20)
      pylab.xlabel("Number of pixels")
      pylab.ylabel("# reflections")
      pylab.savefig(join(self.directory, "n_foreground_hist.png"))
      pylab.clf()


class ReferenceProfileAnalyser(object):
  ''' Analyse the reference profiles. '''

  def __init__(self, directory):
    ''' Set up the directory. '''
    from os.path import join

    # Set the directory
    self.directory = join(directory, "reference")
    ensure_directory(self.directory)

    # Set the required fields
    self.required = [
      "intensity.prf.value",
      "intensity.prf.variance",
      "xyzcal.px",
      "profile.correlation",
    ]

  def __call__(self, rlist):
    ''' Analyse the reference profiles. '''
    from dials.util.command_line import Command

    # Check we have the required fields
    print "Analysing reference profiles"
    if not ensure_required(rlist, self.required):
      return

    # Select only integrated reflections
    Command.start(" Selecting only integated reflections")
    mask = rlist.get_flags(rlist.flags.integrated)
    if mask.count(True) == 0:
      return

    rlist = rlist.select(mask)
    Command.end(" Selected %d integrated reflections" % len(rlist))

    # Analyse distribution of reference spots
    print " Analysing reference profile distribution vs x/y"
    self.reference_xy(rlist)
    print " Analysing reference profile distribution vs z"
    self.reference_z(rlist)

    # Look at correlations between profiles
    def ideal_correlations(filename, rlist):
      ''' Call for reference spots and all reflections. '''

      print " Analysing reflection profile correlations"
      self.ideal_reflection_corr_hist(rlist, filename)

      print " Analysing reflection profile correlations vs x/y"
      self.ideal_reflection_corr_vs_xy(rlist, filename)

      print " Analysing reflection profile correlations vs z"
      self.ideal_reflection_corr_vs_z(rlist, filename)

      print " Analysing reflection profile correlations vs I/Sigma"
      self.ideal_reflection_corr_vs_ios(rlist, filename)

    # Look at correlations between profiles
    def correlations(filename, rlist):
      ''' Call for reference spots and all reflections. '''

      print " Analysing reflection profile correlations"
      self.reflection_corr_hist(rlist, filename)

      print " Analysing reflection profile correlations vs x/y"
      self.reflection_corr_vs_xy(rlist, filename)

      print " Analysing reflection profile correlations vs z"
      self.reflection_corr_vs_z(rlist, filename)

      print " Analysing reflection profile correlations vs I/Sigma"
      self.reflection_corr_vs_ios(rlist, filename)

    mask = rlist.get_flags(rlist.flags.reference_spot)
    correlations("reference",  rlist.select(mask))
    correlations("reflection", rlist)
    ideal_correlations("reference", rlist.select(mask))
    ideal_correlations("reflection", rlist)

  def reference_xy(self, rlist):
    ''' Analyse the distribution of reference profiles. '''
    from os.path import join
    mask = rlist.get_flags(rlist.flags.reference_spot)
    rlist = rlist.select(mask)
    x, y, z = rlist['xyzcal.px'].parts()
    pylab.title("Reference profiles binned in X/Y")
    cax = pylab.hexbin(x, y, gridsize=100)
    pylab.xlabel("x")
    pylab.ylabel("y")
    cbar = pylab.colorbar(cax)
    cbar.ax.set_ylabel("# reflections")
    pylab.savefig(join(self.directory, "reference_xy.png"))
    pylab.clf()

  def reference_z(self, rlist):
    ''' Analyse the distribution of reference profiles. '''
    from os.path import join
    corr = rlist['profile.correlation']
    x, y, z = rlist['xyzcal.px'].parts()
    pylab.title("Reference profiles binned in Z")
    cax = pylab.hist(z, bins=20)
    pylab.xlabel("z")
    pylab.ylabel("# reflections")
    pylab.savefig(join(self.directory, "reference_z.png"))
    pylab.clf()

  def reflection_corr_hist(self, rlist, filename):
    ''' Analyse the correlations. '''
    from os.path import join
    corr = rlist['profile.correlation']
    pylab.title("Reflection correlations histogram")
    pylab.hist(corr, bins=20)
    pylab.xlabel("Correlation with reference profile")
    pylab.ylabel("# reflections")
    pylab.savefig(join(self.directory, "%s_corr_hist" % filename))
    pylab.clf()

  def reflection_corr_vs_xy(self, rlist, filename):
    ''' Analyse the correlations. '''
    from os.path import join
    corr = rlist['profile.correlation']
    x, y, z = rlist['xyzcal.px'].parts()
    pylab.title("Reflection correlations binned in X/Y")
    cax = pylab.hexbin(x, y, C=corr, gridsize=100, vmin=0.0, vmax=1.0)
    cbar = pylab.colorbar(cax)
    pylab.xlabel("x")
    pylab.ylabel("y")
    cbar.ax.set_ylabel("Correlation with reference profile")
    pylab.savefig(join(self.directory, "%s_corr_vs_xy.png" % filename))
    pylab.clf()

  def reflection_corr_vs_z(self, rlist, filename):
    ''' Analyse the correlations. '''
    from os.path import join
    corr = rlist['profile.correlation']
    x, y, z = rlist['xyzcal.px'].parts()
    pylab.title("Reflection correlations vs Z")
    cax = pylab.hexbin(z, corr, gridsize=100)
    cbar = pylab.colorbar(cax)
    pylab.xlabel("z")
    pylab.ylabel("Correlation with reference profile")
    cbar.ax.set_ylabel("# reflections")
    pylab.savefig(join(self.directory, "%s_corr_vs_z.png" % filename))
    pylab.clf()

  def reflection_corr_vs_ios(self, rlist, filename):
    ''' Analyse the correlations. '''
    from dials.array_family import flex
    from os.path import join
    corr = rlist['profile.correlation']
    I = rlist['intensity.prf.value']
    I_sig = flex.sqrt(rlist['intensity.prf.variance'])
    mask = I_sig > 0
    I = I.select(mask)
    I_sig = I_sig.select(mask)
    corr = corr.select(mask)
    I_over_S = I / I_sig
    mask = I_over_S > 0.1
    I_over_S = I_over_S.select(mask)
    corr = corr.select(mask)
    pylab.title("Reflection correlations vs Log I/Sigma")
    cax = pylab.hexbin(flex.log(I_over_S), corr, gridsize=100)
    cbar = pylab.colorbar(cax)
    pylab.xlabel("Log I/Sigma")
    pylab.ylabel("Correlation with reference profile")
    cbar.ax.set_ylabel("# reflections")
    pylab.savefig(join(self.directory, "%s_corr_vs_ios.png" % filename))
    pylab.clf()

  def ideal_reflection_corr_hist(self, rlist, filename):
    ''' Analyse the correlations. '''
    from os.path import join
    if 'correlation.ideal.profile' in rlist:
      corr = rlist['correlation.ideal.profile']
      pylab.title("Reflection correlations histogram")
      pylab.hist(corr, bins=20)
      pylab.xlabel("Correlation with reference profile")
      pylab.ylabel("# reflections")
      pylab.savefig(join(self.directory, "ideal_%s_corr_hist" % filename))
      pylab.clf()

  def ideal_reflection_corr_vs_xy(self, rlist, filename):
    ''' Analyse the correlations. '''
    from os.path import join
    if 'correlation.ideal.profile' in rlist:
      corr = rlist['correlation.ideal.profile']
      x, y, z = rlist['xyzcal.px'].parts()
      pylab.title("Reflection correlations binned in X/Y")
      cax = pylab.hexbin(x, y, C=corr, gridsize=100, vmin=0.0, vmax=1.0)
      cbar = pylab.colorbar(cax)
      pylab.xlabel("x")
      pylab.ylabel("y")
      cbar.ax.set_ylabel("Correlation with reference profile")
      pylab.savefig(join(self.directory, "ideal_%s_corr_vs_xy.png" % filename))
      pylab.clf()

  def ideal_reflection_corr_vs_z(self, rlist, filename):
    ''' Analyse the correlations. '''
    from os.path import join
    if 'correlation.ideal.profile' in rlist:
      corr = rlist['correlation.ideal.profile']
      x, y, z = rlist['xyzcal.px'].parts()
      pylab.title("Reflection correlations vs Z")
      cax = pylab.hexbin(z, corr, gridsize=100)
      cbar = pylab.colorbar(cax)
      pylab.xlabel("z")
      pylab.ylabel("Correlation with reference profile")
      cbar.ax.set_ylabel("# reflections")
      pylab.savefig(join(self.directory, "ideal_%s_corr_vs_z.png" % filename))
      pylab.clf()

  def ideal_reflection_corr_vs_ios(self, rlist, filename):
    ''' Analyse the correlations. '''
    from dials.array_family import flex
    from os.path import join
    if 'correlation.ideal.profile' in rlist:
      corr = rlist['correlation.ideal.profile']
      I = rlist['intensity.prf.value']
      I_sig = flex.sqrt(rlist['intensity.prf.variance'])
      mask = I_sig > 0
      I = I.select(mask)
      I_sig = I_sig.select(mask)
      corr = corr.select(mask)
      I_over_S = I / I_sig
      mask = I_over_S > 0.1
      I_over_S = I_over_S.select(mask)
      corr = corr.select(mask)
      pylab.title("Reflection correlations vs Log I/Sigma")
      cax = pylab.hexbin(flex.log(I_over_S), corr, gridsize=100)
      cbar = pylab.colorbar(cax)
      pylab.xlabel("Log I/Sigma")
      pylab.ylabel("Correlation with reference profile")
      cbar.ax.set_ylabel("# reflections")
      pylab.savefig(join(self.directory, "ideal_%s_corr_vs_ios.png" % filename))
      pylab.clf()


class Analyser(object):
  ''' Helper class to do all the analysis. '''

  def __init__(self, directory):
    ''' Setup the analysers. '''
    from os.path import join
    directory = join(directory, "analysis")
    self.analysers = [
      CentroidAnalyser(directory),
      BackgroundAnalyser(directory),
      IntensityAnalyser(directory),
      ReferenceProfileAnalyser(directory),
    ]

  def __call__(self, rlist):
    ''' Do all the analysis. '''
    from copy import deepcopy
    for analyse in self.analysers:
      analyse(deepcopy(rlist))


class Script(object):
  ''' A class to encapsulate the script. '''

  def __init__(self):
    ''' Initialise the script. '''
    from dials.util.options import OptionParser
    from libtbx.phil import parse
    import libtbx.load_env

    # Create the phil parameters
    phil_scope = parse('''
      output {
        directory = .
          .type = str
          .help = "The directory to store the results"
      }
    ''')

    # Create the parser
    usage = "usage: %s [options] reflections.pickle" % libtbx.env.dispatcher_name
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      read_reflections=True)

    self.parser.add_option(
      '--xkcd',
      action='store_true',
      dest='xkcd',
      default=False,
      help='Special drawing mode')

  def run(self):
    ''' Run the script. '''
    from dials.array_family import flex
    from dials.util.command_line import Command

    # Parse the command line arguments
    params, options = self.parser.parse_args(show_diff_phil=True)

    if options.xkcd:
      from matplotlib import pyplot
      pyplot.xkcd()

    # Shoe the help
    if len(params.input.reflections) != 1:
      self.parser.print_help()
      exit(0)

    # Analyse the reflections
    analyse = Analyser(params.output.directory)
    analyse(params.input.reflections[0].data)


if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
