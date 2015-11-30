# -*- coding: utf-8 -*-
#!/usr/bin/env python
#
# analyse_output.py
#
#  Copyright (C) 2015 Diamond Light Source
#
#  Author: James Parkhurst, Richard Gildea
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from __future__ import division
import copy
import math

from libtbx.containers import OrderedDict
from dials.array_family import flex

RAD2DEG = 180/math.pi

help_message = '''

Generate a number of analysis plots from input integrated or indexed reflections.

Examples::

  dials.analyse_output indexed.pickle

  dials.analyse_output refined.pickle

  dials.analyse_output integrated.pickle

'''

import libtbx.phil

# Create the phil parameters
phil_scope = libtbx.phil.parse('''
  output {
    html = dials-report.html
      .type = path
      .help = "The name of the output html file"
    json = None
      .type = path
      .help = "The name of the optional json file containing the plot data"
  }
  grid_size = Auto
    .type = ints(size=2)
  pixels_per_bin = 40
    .type = int(value_min=1)

  centroid_diff_max = None
    .help = "Magnitude in pixels of shifts mapped to the extreme colours"
            "in the heatmap plots centroid_diff_x and centroid_diff_y"
    .type = float
    .expert_level = 1

  orientation_decomposition
    .help = "Options determining how the orientation matrix"
            "decomposition is done. The axes about which to decompose"
            "the matrix into three rotations are chosen here, as well"
            "as whether the rotations are relative to the reference"
            "orientation, taken from the static crystal model"
  {
    e1 = 1. 0. 0.
      .type = floats(size = 3)

    e2 = 0. 1. 0.
      .type = floats(size = 3)

    e3 = 0. 0. 1.
      .type = floats(size = 3)

    relative_to_static_orientation = True
      .type = bool
  }
''')


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

def determine_grid_size(rlist, grid_size=None):
  from libtbx import Auto
  panel_ids = rlist['panel']
  n_panels = flex.max(panel_ids) + 1
  if grid_size is not None and grid_size is not Auto:
    assert (grid_size[0] * grid_size[1]) >= n_panels, (n_panels)
    return grid_size
  n_cols = int(math.floor(math.sqrt(n_panels)))
  n_rows = int(math.ceil(n_panels / n_cols))
  return n_cols, n_rows


class per_panel_plot(object):

  title = None
  filename = None
  cbar_ylabel = None
  xlabel = 'x'
  ylabel = 'y'

  def __init__(self, rlist, directory, grid_size=None, pixels_per_bin=10):

    from os.path import join

    min_x, max_x, min_y, max_y = self.get_min_max_xy(rlist)
    panel_ids = rlist['panel']
    crystal_ids = rlist['id']
    n_crystals = flex.max(crystal_ids) + 1
    n_panels = flex.max(panel_ids) + 1

    n_cols, n_rows = determine_grid_size(
      rlist, grid_size=grid_size)

    from matplotlib import pyplot

    for i_crystal in range(n_crystals):
      if n_crystals > 1:
        suffix = '_%i' %i_crystal
      else:
        suffix = ''
      crystal_sel = (crystal_ids == i_crystal)
      fig, axes = pyplot.subplots(
        n_rows, n_cols)

      if n_panels == 1:
        axes = [[axes]]
      elif n_cols == 1:
        axes = [[ax] for ax in axes]
      elif n_rows == 1:
        axes_x = [axes]

      self.gridsize = tuple(
        int(math.ceil(i))
        for i in (max_x/pixels_per_bin, max_y/pixels_per_bin))

      clim = (1e8, 1e-8)

      plots = []

      i_panel = 0
      for i_row in range(n_rows):
        for i_col in range(n_cols):

          panel_sel = (panel_ids == i_panel)
          sel = (panel_sel & crystal_sel)
          i_panel += 1

          if n_panels > 1:
            axes[i_row][i_col].set_title('Panel %d' %i_panel)
            axes[i_row][i_col].set_title('Panel %d' %i_panel)

          if (i_row+1) == n_rows:
            axes[i_row][i_col].set_xlabel(self.xlabel)
          else:
            pyplot.setp(axes[i_row][i_col].get_xticklabels(), visible=False)

          if i_col == 0:
            axes[i_row][i_col].set_ylabel(self.ylabel)
          else:
            pyplot.setp(axes[i_row][i_col].get_yticklabels(), visible=False)

          if sel.count(True) > 0:
            rlist_sel = rlist.select(sel)
            if len(rlist_sel) <= 1:
              ax = pyplot.scatter([],[]) # create empty plot
            else:
              ax = self.plot_one_panel(axes[i_row][i_col], rlist_sel)
              clim = (min(clim[0], ax.get_clim()[0]),
                      max(clim[1], ax.get_clim()[1]))
            plots.append(ax)

          axes[i_row][i_col].set_xlim(min_x, max_x)
          axes[i_row][i_col].set_ylim(min_y, max_y)
          axes[i_row][i_col].axes.set_aspect('equal')
          axes[i_row][i_col].invert_yaxis()

      for p in plots:
        p.set_clim(clim)

      default_size = fig.get_size_inches()
      fig.set_size_inches((n_cols*default_size[0], n_rows*default_size[1]))

      #pyplot.tight_layout()
      if self.cbar_ylabel is not None:
        cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
        cbar = fig.colorbar(ax, cax=cax)
        cbar.ax.set_ylabel(self.cbar_ylabel, fontsize=n_cols*10)
        cbar.ax.tick_params(labelsize=n_cols*8)
        if n_panels > 1:
          fig.subplots_adjust(hspace=0.1/n_rows, right=0.8)

      if self.title is not None:
        fig.suptitle(self.title, fontsize=n_cols*12)
      fig.savefig(join(directory, self.filename))
      fig.set_size_inches(default_size)
      pyplot.close()

  def get_min_max_xy(self, rlist):
    xc, yc, zc = rlist['xyzcal.px'].parts()
    xo, yo, zo = rlist['xyzobs.px.value'].parts()

    min_x = math.floor(min(flex.min(xc), flex.min(xo)))
    min_y = math.floor(min(flex.min(yc), flex.min(yo)))
    max_x = math.ceil(max(flex.max(xc), flex.max(xo)))
    max_y = math.ceil(max(flex.max(yc), flex.max(yo)))
    return min_x, max_x, min_y, max_y

  def plot_one_panel(self, ax, rlist):
    raise NotImplementedError()


class ScanVaryingCrystalAnalyser(object):
  ''' Analyse a scan-varying crystal. '''

  def __init__(self, orientation_decomposition):
    # Decomposition axes
    self._e1 = orientation_decomposition.e1
    self._e2 = orientation_decomposition.e2
    self._e3 = orientation_decomposition.e3
    self._relative_to_static_orientation \
      = orientation_decomposition.relative_to_static_orientation
    self._debug = False

  def __call__(self, experiments):
    ''' Analyse the strong spots. '''
    from dials.util.command_line import Command

    # Check we have the required fields
    print "Analysing scan-varying crystal model"

    d = OrderedDict()

    if experiments is not None and len(experiments):
      d.update(self.plot_cell(experiments))
      d.update(self.plot_orientation(experiments))
    return {'scan_varying': d}

  def plot_cell(self, experiments):
    ''' Analyse the scan-varying cell parameters. '''

    # cell plot
    dat = []
    for iexp, exp in enumerate(experiments):

      crystal = exp.crystal
      scan = exp.scan

      if crystal.num_scan_points == 0:
        print "Ignoring scan-static crystal"
        continue

      scan_pts = range(crystal.num_scan_points)
      cells = [crystal.get_unit_cell_at_scan_point(t) for t in scan_pts]
      cell_params = [e.parameters() for e in cells]
      a,b,c,aa,bb,cc = zip(*cell_params)
      aa = list(round(i, ndigits=6) for i in aa)
      bb = list(round(i, ndigits=6) for i in bb)
      cc = list(round(i, ndigits=6) for i in cc)
      phi = [scan.get_angle_from_array_index(t) for t in scan_pts]
      vol = [e.volume() for e in cells]
      cell_dat = {'phi':phi,
                  'a':a,
                  'b':b,
                  'c':c,
                  'alpha':aa,
                  'beta':bb,
                  'gamma':cc,
                  'volume':vol}
      if self._debug:
        print "Crystal in Experiment {0}".format(iexp)
        print "Phi\ta\tb\tc\talpha\tbeta\tgamma\tVolume"
        msg = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}"
        line_dat = zip(phi, a, b, c, aa, bb, cc, vol)
        for line in line_dat:
          print msg.format(*line)
      dat.append(cell_dat)

    d = {
      'scan_varying_cell': {
        'data': [],
        'layout': {
          'title': 'Scan-varying cell parameters',
          'xaxis3': {
            'domain': [0, 1],
            'anchor': 'y7',
            'title': 'rotation angle (°)',
          },
          'xaxis2': {
            'domain': [0.55, 1],
            'anchor': 'y6',
          },
          'xaxis': {
            'domain': [0, 0.45],
            'anchor': 'y3',
          },

          'yaxis7': {
            'domain': [0.0, 0.2],
            'anchor': 'x3',
            'nticks': 3,
          },
          'yaxis6': {
            'domain': [0.3, 0.5],
            'anchor': 'x2',
            'nticks': 3,
          },
          'yaxis5': {
            'domain': [0.55, 0.75],
            'anchor': 'x2',
            'nticks': 3,
          },
          'yaxis4': {
            'domain': [0.8, 1],
            'anchor': 'x2',
            'nticks': 3,
          },
          'yaxis3': {
            'domain': [0.3, 0.5],
            'anchor': 'x',
            'nticks': 3,
          },
          'yaxis2': {
            'domain': [0.55, 0.75],
            'anchor': 'x',
            'nticks': 3,
          },
          'yaxis': {
            'domain': [0.8, 1],
            'nticks': 3,
          },

        },
      },
    }

    for cell_dat in dat:
      d['scan_varying_cell']['data'].extend([
        {
          'x': cell_dat['phi'],
          'y': cell_dat['a'],
          'type': 'scatter',
          'name': 'a (Å)',
        },
        {
          'x': cell_dat['phi'],
          'y': cell_dat['b'],
          'type': 'scatter',
          'name': 'b (Å)',
          'xaxis': 'x',
          'yaxis': 'y2',
        },
        {
          'x': cell_dat['phi'],
          'y': cell_dat['c'],
          'type': 'scatter',
          'name': 'c (Å)',
          'xaxis': 'x',
          'yaxis': 'y3',
        },
        {
          'x': cell_dat['phi'],
          'y': cell_dat['alpha'],
          'type': 'scatter',
          'name': 'α (°)',
          'xaxis': 'x2',
          'yaxis': 'y4',
        },
        {
          'x': cell_dat['phi'],
          'y': cell_dat['beta'],
          'type': 'scatter',
          'name': 'β (°)',
          'xaxis': 'x2',
          'yaxis': 'y5',
        },
        {
          'x': cell_dat['phi'],
          'y': cell_dat['gamma'],
          'type': 'scatter',
          'name': 'γ (°)',
          'xaxis': 'x2',
          'yaxis': 'y6',
        },
        {
          'x': cell_dat['phi'],
          'y': cell_dat['volume'],
          'type': 'scatter',
          'name': 'volume (Å^3)',
          'xaxis': 'x3',
          'yaxis': 'y7',
        }])
    if not len(dat):
      return {}
    return d

  def plot_orientation(self, experiments):
    from dials.algorithms.refinement.rotation_decomposition import \
      solve_r3_rotation_for_angles_given_axes

    # orientation plot
    dat = []
    for iexp, exp in enumerate(experiments):

      crystal = exp.crystal
      scan = exp.scan

      if crystal.num_scan_points == 0:
        print "Ignoring scan-static crystal"
        continue

      scan_pts = range(crystal.num_scan_points)
      phi = [scan.get_angle_from_array_index(t) for t in scan_pts]
      Umats = [crystal.get_U_at_scan_point(t) for t in scan_pts]
      if self._relative_to_static_orientation:
        # factor out static U
        Uinv = crystal.get_U().inverse()
        Umats = [U*Uinv for U in Umats]
      # NB e3 and e1 definitions for the crystal are swapped compared
      # with those used inside the solve_r3_rotation_for_angles_given_axes
      # method
      angles = [solve_r3_rotation_for_angles_given_axes(U,
        self._e3, self._e2, self._e1, deg=True) for U in Umats]
      phi3, phi2, phi1 = zip(*angles)
      angle_dat = {'phi':phi,
                   'phi3':phi3,
                   'phi2':phi2,
                   'phi1':phi1}
      if self._debug:
        print "Crystal in Experiment {0}".format(iexp)
        print "Image\tphi3\tphi2\tphi1"
        msg = "{0}\t{1}\t{2}\t{3}"
        line_dat = zip(phi, phi3, phi2, phi1)
        for line in line_dat:
          print msg.format(*line)
      dat.append(angle_dat)

    d = {
      'scan_varying_orientation': {
        'data': [],
        'layout': {
          'title': 'Scan-varying orientation parameters',
          'xaxis3': {
            'domain': [0, 1],
            'anchor': 'y3',
            'title': 'rotation angle (°)',
          },
          'xaxis2': {
            'domain': [0, 1],
            'anchor': 'y3',
          },
          'xaxis': {
            'domain': [0, 1],
            'anchor': 'y3',
          },

          'yaxis3': {
            'domain': [0, 0.3],
            'anchor': 'x3',
          },
          'yaxis2': {
            'domain': [0.35, 0.65],
            'anchor': 'x2',
          },
          'yaxis': {
            'domain': [0.7, 1],
          },

        },
      },
    }

    for ori in dat:
      d['scan_varying_orientation']['data'].extend([
        {
          'x': ori['phi'],
          'y': ori['phi1'],
          'type': 'scatter',
          'name': 'Φ1 (°)',
        },
        {
          'x': ori['phi'],
          'y': ori['phi2'],
          'type': 'scatter',
          'name': 'Φ2 (°)',
          'xaxis': 'x2',
          'yaxis': 'y2',
        },
        {
          'x': ori['phi'],
          'y': ori['phi3'],
          'type': 'scatter',
          'name': 'Φ3 (°)',
          'xaxis': 'x3',
          'yaxis': 'y3',
        }])
    if not len(dat):
      return {}
    return d


class StrongSpotsAnalyser(object):
  ''' Analyse a list of strong spots. '''

  def __init__(self):
    from os.path import join

    # Set the required fields
    self.required = [
      "xyzobs.px.value",
      "panel",
    ]

  def __call__(self, rlist):
    ''' Analyse the strong spots. '''
    from dials.util.command_line import Command

    # Check we have the required fields
    print "Analysing strong spots"
    if not ensure_required(rlist, self.required):
      return {'strong': {}}

    # Remove I_sigma <= 0
    selection = rlist['intensity.sum.variance'] <= 0
    if selection.count(True) > 0:
      rlist.del_selected(selection)
      print ' Removing %d reflections with variance <= 0' % \
        selection.count(True)

    if 'flags' in rlist:
      # Select only strong reflections
      Command.start(" Selecting only strong reflections")
      mask = rlist.get_flags(rlist.flags.strong)
      if mask.count(True) > 0:
        threshold = 10
        rlist = rlist.select(mask)
      Command.end(" Selected %d strong reflections" % len(rlist))

    d = OrderedDict()
    # Look at distribution of spot counts
    d.update(self.spot_count_per_image(rlist))
    self.spot_count_per_panel(rlist)
    return {'strong': d}

  def spot_count_per_image(self, rlist):
    ''' Analyse the spot count per image. '''
    from os.path import join
    x,y,z = rlist['xyzobs.px.value'].parts()
    max_z = int(math.ceil(flex.max(z)))

    ids = rlist['id']
    spot_count_per_image = []
    for j in range(flex.max(ids)+1):
      spot_count_per_image.append([])
      zsel = z.select(ids == j)
      for i in range(max_z):
        sel = (zsel >= i) & (zsel < (i+1))
        spot_count_per_image[j].append(sel.count(True))

    d = {
      'spot_count_per_image': {
        'data': [],
        'layout': {
          'title': 'Spot count per image',
          'xaxis': {'title': 'Image'},
          'yaxis': {
            'title': 'Spot count',
            'rangemode': 'tozero'
          },
        },
      },
    }
    for j in range(flex.max(ids)+1):
      d['spot_count_per_image']['data'].append({
        'x': list(range(len(spot_count_per_image[j]))),
        'y': spot_count_per_image[j],
        'type': 'scatter',
        'name': '#spots',
        'opacity': 0.4,
      })
    return d

  def spot_count_per_panel(self, rlist):
    ''' Analyse the spot count per panel. '''
    from os.path import join
    panel = rlist['panel']
    if flex.max(panel) == 0:
      # only one panel, don't bother generating a plot
      return

    n_panels = int(flex.max(panel))
    spot_count_per_panel = flex.int()
    for i in range(n_panels):
      sel = (panel >= i) & (panel < (i+1))
      spot_count_per_panel.append(sel.count(True))

    from matplotlib import pyplot
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    ax.set_title("Spot count per panel")
    ax.scatter(
      list(range(len(spot_count_per_panel))), spot_count_per_panel,
      s=10, color='blue', marker='o', alpha=0.4)
    ax.set_xlabel("Panel #")
    ax.set_ylabel("# spots")
    pyplot.savefig(join(self.directory, "spots_per_panel.png"))
    pyplot.close()


class CentroidAnalyser(object):
  ''' Analyse the reflection centroids. '''

  def __init__(self, grid_size=None, pixels_per_bin=10,
    centroid_diff_max=1.5):
    from os.path import join

    self.grid_size = grid_size
    self.pixels_per_bin = pixels_per_bin
    self.centroid_diff_max = centroid_diff_max

    # Set the required fields
    self.required = [
      "intensity.sum.value",
      "intensity.sum.variance",
      "xyzcal.px",
      "xyzobs.px.value",
      "xyzcal.mm",
      "xyzobs.mm.value",

    ]

  def __call__(self, rlist):
    ''' Analyse the reflection centroids. '''
    from dials.util.command_line import Command

    # Check we have the required fields
    print "Analysing reflection centroids"
    if not ensure_required(rlist, self.required):
      return {'centroid': {}}

    # Remove I_sigma <= 0
    selection = rlist['intensity.sum.variance'] <= 0
    if selection.count(True) > 0:
      rlist.del_selected(selection)
      print ' Removing %d reflections with variance <= 0' % \
        selection.count(True)

    # Remove partial reflections as their observed centroids won't be accurate
    if 'partiality' in rlist:
      selection = rlist['partiality'] < 0.99
      if selection.count(True) > 0 and selection.count(True) < selection.size():
        rlist.del_selected(selection)
        print ' Removing %d partial reflections' % \
          selection.count(True)

    # Select only integrated reflections
    Command.start(" Selecting only summation-integated reflections")
    mask = rlist.get_flags(rlist.flags.integrated_sum)
    if mask.count(True) > 0:
      threshold = 10
      rlist = rlist.select(mask)
      Command.end(" Selected %d summation-integrated reflections" % len(rlist))
    else:
      # Select only those reflections used in refinement
      threshold = 0
      mask = rlist.get_flags(rlist.flags.used_in_refinement)
      rlist = rlist.select(mask)
      Command.end(" Selected %d refined reflections" % len(rlist))

    d = OrderedDict()

    # Look at differences in calculated/observed position
    print " Analysing centroid differences with I/Sigma > %s" %threshold
    d.update(self.centroid_diff_hist(rlist, threshold))
    print " Analysing centroid differences in x/y with I/Sigma > %s" %threshold
    d.update(self.centroid_diff_xy(rlist, threshold))
    d.update(self.centroid_xy_xz_zy_residuals(rlist, threshold))
    print " Analysing centroid differences in z with I/Sigma > %s" %threshold
    d.update(self.centroid_diff_z(rlist, threshold))
    print " Analysing centroid differences vs phi with I/Sigma > %s" %threshold
    d.update(self.centroid_mean_diff_vs_phi(rlist, threshold))
    return {'centroid': d}

  def centroid_diff_hist(self, rlist, threshold):
    ''' Analyse the correlations. '''
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
    hist = flex.histogram(diff, n_slots=20)

    d = {
      'centroid_difference_histogram': {
        'data': [{
          'x': list(hist.slot_centers()),
          'y': list(hist.slots()),
          'type': 'bar',
          'name': 'centroid_differences',
        }],
        'layout': {
          'title': 'Difference between observed and calculated centroids',
          'xaxis': {'title': 'Difference in position'},
          'yaxis': {
            'title': 'Number of reflections',
          },
          'bargap': 0,
        },
      },
    }
    return d

  def centroid_diff_xy(self, rlist, threshold):
    ''' Look at the centroid difference in x, y '''
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

    d = OrderedDict()

    nbinsx, nbinsy = tuple(
        int(math.ceil(i))
        for i in (flex.max(xc)/self.pixels_per_bin,
                  flex.max(yc)/self.pixels_per_bin))

    xc = xc.as_numpy_array()
    yc = yc.as_numpy_array()
    xd = xd.as_numpy_array()
    yd = yd.as_numpy_array()

    import numpy as np
    H, xedges, yedges = np.histogram2d(xc, yc, bins=(nbinsx, nbinsy))
    H1, xedges, yedges = np.histogram2d(xc, yc, bins=(nbinsx, nbinsy), weights=xd)
    H2, xedges, yedges = np.histogram2d(xc, yc, bins=(nbinsx, nbinsy), weights=yd)

    nonzeros = np.nonzero(H)
    z1 = np.empty(H.shape)
    z1[:] = np.NAN
    z1[nonzeros] = (H1[nonzeros]/H[nonzeros])
    z2 = np.empty(H.shape)
    z2[:] = np.NAN
    z2[nonzeros] = (H2[nonzeros]/H[nonzeros])

    d.update({
      'centroid_differences_x': {
        'data': [{
          'name': 'centroid_differences_x',
          'x': xedges.tolist(),
          'y': xedges.tolist(),
          'z': z1.transpose().tolist(),
          'type': 'heatmap',
          'colorbar': {
            'title': 'Difference in X position',
            'titleside': 'right',
          },
          'colorscale': 'Jet',
        }],
        'layout': {
          'title': 'Difference between observed and calculated centroids in X',
          'xaxis': {
            'domain': [0, 0.85],
            'title': 'X',
            'showgrid': False,
          },
          'yaxis': {
            'title': 'Y',
            'autorange': 'reversed',
            'showgrid': False,
          },
          'width': 500,
          'height': 450,
        },
      },
    })

    d.update({
      'centroid_differences_y': {
        'data': [{
          'name': 'centroid_differences_y',
          'x': xedges.tolist(),
          'y': xedges.tolist(),
          'z': z2.transpose().tolist(),
          'type': 'heatmap',
          'colorbar': {
            'title': 'Difference in Y position',
            'titleside': 'right',
          },
          'colorscale': 'Jet',
        }],
        'layout': {
          'title': 'Difference between observed and calculated centroids in Y',
          'xaxis': {
            'domain': [0, 0.85],
            'title': 'X',
            'showgrid': False,
          },
          'yaxis': {
            'title': 'Y',
            'autorange': 'reversed',
            'showgrid': False,
          },
          'width': 500,
          'height': 450,
        },
      },
    })
    return d

  def centroid_diff_z(self, rlist, threshold):
    ''' Look at the centroid difference in z '''
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

    import numpy as np
    H, xedges, yedges = np.histogram2d(zc, zd, bins=(100, 100))

    return {
      'centroid_differences_z': {
        'data': [{
          'x': xedges.tolist(),
          'y': xedges.tolist(),
          'z': H.transpose().tolist(),
          'type': 'heatmap',
          'name': 'centroid_differences_z',
          'colorbar': {
            'title': 'Number of reflections',
            'titleside': 'right',
          },
          'colorscale': 'Jet',
        }],
        'layout': {
          'title': 'Difference between observed and calculated centroids in Z',
          'xaxis': {
            'title': 'Z',
            'showgrid': False,
          },
          'yaxis': {
            'title': 'Difference in Z position',
            'showgrid': False,
          },
        },
      },
    }

  def centroid_mean_diff_vs_phi(self, rlist, threshold):
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
    dphi = (zc - zo) * RAD2DEG

    mean_residuals_x = flex.double()
    mean_residuals_y = flex.double()
    mean_residuals_phi = flex.double()
    rmsd_x = flex.double()
    rmsd_y = flex.double()
    rmsd_phi = flex.double()
    frame = []
    phi_obs_deg = RAD2DEG * zo
    phi = []

    for i_phi in range(int(math.floor(flex.min(phi_obs_deg))),
                   int(math.ceil(flex.max(phi_obs_deg)))):
      sel = (phi_obs_deg >= i_phi) & (phi_obs_deg < (i_phi+1))
      if sel.count(True) == 0:
        continue
      mean_residuals_x.append(flex.mean(dx.select(sel)))
      mean_residuals_y.append(flex.mean(dy.select(sel)))
      mean_residuals_phi.append(flex.mean(dphi.select(sel)))
      rmsd_x.append(math.sqrt(flex.mean_sq(dx.select(sel))))
      rmsd_y.append(math.sqrt(flex.mean_sq(dy.select(sel))))
      rmsd_phi.append(math.sqrt(flex.mean_sq(dphi.select(sel))))
      phi.append(i_phi)

    d = {
      'centroid_mean_differences_vs_phi': {
        'data': [
          {
            'x': list(phi),
            'y': list(mean_residuals_x),
            'type': 'scatter',
            'name': 'mean_dx',
          },
          {
            'x': list(phi),
            'y': list(mean_residuals_y),
            'type': 'scatter',
            'name': 'mean_dy',
            'xaxis': 'x2',
            'yaxis': 'y2',
          },
          {
            'x': list(phi),
            'y': list(mean_residuals_phi),
            'type': 'scatter',
            'name': 'mean_dphi',
            'xaxis': 'x3',
            'yaxis': 'y3',
          },
        ],
        'layout': {
          'title': 'Difference between observed and calculated centroids vs z',
          'yaxis3': {'domain': [0, 0.266]},
          #'legend': {'traceorder': 'reversed'},
          'xaxis3': {'anchor': 'y3'},
          'xaxis2': {'anchor': 'y2'},
          'yaxis2': {'domain': [0.366, 0.633]},
          'yaxis': {'domain': [0.733, 1]}
        },
      },
      'centroid_rmsd_vs_phi': {
        'data': [
          {
            'x': list(phi),
            'y': list(rmsd_x),
            'type': 'scatter',
            'name': 'rmsd_dx',
          },
          {
            'x': list(phi),
            'y': list(rmsd_y),
            'type': 'scatter',
            'name': 'rmsd_dy',
            'xaxis': 'x2',
            'yaxis': 'y2',
          },
          {
            'x': list(phi),
            'y': list(rmsd_phi),
            'type': 'scatter',
            'name': 'rmsd_dphi',
            'xaxis': 'x3',
            'yaxis': 'y3',
          },
        ],
        'layout': {
          'title': 'RMSD between observed and calculated centroids vs phi',
          'yaxis3': {'domain': [0, 0.266]},
          #'legend': {'traceorder': 'reversed'},
          'xaxis3': {'anchor': 'y3'},
          'xaxis2': {'anchor': 'y2'},
          'yaxis2': {'domain': [0.366, 0.633]},
          'yaxis': {'domain': [0.733, 1]}
        },
      },
    }
    return d

  def centroid_xy_xz_zy_residuals(self, rlist, threshold):
    from os.path import join
    import math
    I = rlist['intensity.sum.value']
    I_sig = flex.sqrt(rlist['intensity.sum.variance'])
    I_over_S = I / I_sig
    mask = I_over_S > threshold
    rlist = rlist.select(mask)
    assert(len(rlist) > 0)

    xc, yc, zc = rlist['xyzcal.px'].parts()
    xo, yo, zo = rlist['xyzobs.px.value'].parts()
    dx = xc - xo
    dy = yc - yo
    dz = zc - zo

    d = OrderedDict()

    import numpy as np
    Hxy, xedges, yedges = np.histogram2d(dx, dy, bins=(50, 50))
    Hzy, zedges, yedges = np.histogram2d(dz, dy, bins=(50, 50))
    Hxz, xedges, zedges = np.histogram2d(dx, dz, bins=(50, 50))

    density_hist_layout = {

      'showlegend': False,
      'autosize': False,
      'width': 600,
      'height': 550,
      'margin': {'t': 50},
      'hovermode': 'closest',
      'bargap': 0,
      'xaxis': {
        'domain': [0, 0.85],
        'showgrid': False,
        'zeroline': True,
        'zerolinewidth': 2,
        'zerolinecolor': '#969696',
        'title': 'X',
      },
      'yaxis': {
        'domain': [0, 0.85],
        'showgrid': False,
        'zeroline': True,
        'zerolinewidth': 2,
        'zerolinecolor': '#969696',
        'title': 'Y',
      },
      'xaxis2': {
        'domain': [0.85, 1],
        'showgrid': False,
        'zeroline': False,
        'bargap': 0,
      },
      'yaxis2': {
        'domain': [0.85, 1],
        'showgrid': False,
        'zeroline': False,
        'bargap': 0,
      }
    }

    histx = flex.histogram(dx, n_slots=100)
    histy = flex.histogram(dy, n_slots=100)
    histz = flex.histogram(dz, n_slots=100)

    d.update({
      'residuals_xy': {
        'data': [
          #{
            #'x': list(dx),
            #'y': list(dy),
            #'mode': 'markers',
            #'name': 'points',
            #'marker': {
              #'color': 'rgb(102,0,0)',
              #'size': 2,
              #'opacity': 0.4
            #},
            #'type': 'scatter',
          #},
          {
            'x': xedges.tolist(),
            'y': yedges.tolist(),
            'z': Hxy.transpose().tolist(),
            'name': 'density',
            #'ncontours': 20,
            'colorscale': 'Hot',
            'reversescale': True,
            'showscale': False,
            'type': 'contour',
            'zsmooth': 'best',
          },
          {
            'x': list(histx.slot_centers()),
            'y': list(histx.slots()),
            'name': 'dx histogram',
            'marker': {'color': 'rgb(102,0,0)'},
            'yaxis': 'y2',
            'type': 'bar',
          },
          {
            'y': list(histy.slot_centers()),
            'x': list(histy.slots()),
            'name': 'dy histogram',
            'marker': {'color': 'rgb(102,0,0)'},
            'xaxis': 'x2',
            'type': 'bar',
            'orientation': 'h',
          },
        ],
        'layout': copy.deepcopy(density_hist_layout),
      }
    })
    d['residuals_xy']['layout']['title'] = 'Centroid residuals in X and Y'

    d.update({
      'residuals_zy': {
        'data': [
          #{
            #'x': list(dz),
            #'y': list(dy),
            #'mode': 'markers',
            #'name': 'points',
            #'marker': {
              #'color': 'rgb(102,0,0)',
              #'size': 2,
              #'opacity': 0.4
            #},
            #'type': 'scatter',
          #},
          {
            'x': zedges.tolist(),
            'y': yedges.tolist(),
            'z': Hzy.transpose().tolist(),
            'name': 'density',
            #'ncontours': 20,
            'colorscale': 'Hot',
            'reversescale': True,
            'showscale': False,
            'type': 'contour',
            'zsmooth': 'best',
          },
          {
            'x': list(histz.slot_centers()),
            'y': list(histz.slots()),
            'name': 'dz histogram',
            'marker': {'color': 'rgb(102,0,0)'},
            'yaxis': 'y2',
            'type': 'bar',
          },
          {
            'y': list(histy.slot_centers()),
            'x': list(histy.slots()),
            'name': 'dy histogram',
            'marker': {'color': 'rgb(102,0,0)'},
            'xaxis': 'x2',
            'type': 'bar',
            'orientation': 'h',
          },
        ],
        'layout': copy.deepcopy(density_hist_layout),
      }
    })
    d['residuals_zy']['layout']['title'] = 'Centroid residuals in Z and Y'
    d['residuals_zy']['layout']['xaxis']['title'] = 'Z'

    d.update({
      'residuals_xz': {
        'data': [
          #{
            #'x': list(dx),
            #'y': list(dz),
            #'mode': 'markers',
            #'name': 'points',
            #'marker': {
              #'color': 'rgb(102,0,0)',
              #'size': 2,
              #'opacity': 0.4
            #},
            #'type': 'scatter',
          #},
          {
            'x': xedges.tolist(),
            'y': zedges.tolist(),
            'z': Hxz.transpose().tolist(),
            'name': 'density',
            #'ncontours': 20,
            'colorscale': 'Hot',
            'reversescale': True,
            'showscale': False,
            'type': 'contour',
            'zsmooth': 'best',
          },
          {
            'x': list(histx.slot_centers()),
            'y': list(histx.slots()),
            'name': 'dx histogram',
            'marker': {'color': 'rgb(102,0,0)'},
            'yaxis': 'y2',
            'type': 'bar',
          },
          {
            'y': list(histz.slot_centers()),
            'x': list(histz.slots()),
            'name': 'dz histogram',
            'marker': {'color': 'rgb(102,0,0)'},
            'xaxis': 'x2',
            'type': 'bar',
            'orientation': 'h',
          },
        ],
        'layout': copy.deepcopy(density_hist_layout),
      }
    })
    d['residuals_xz']['layout']['title'] = 'Centroid residuals in X and Z'
    d['residuals_xz']['layout']['yaxis']['title'] = 'Z'

    return d


class BackgroundAnalyser(object):
  ''' Analyse the background. '''

  def __init__(self, grid_size=None, pixels_per_bin=10):
    from os.path import join

    self.grid_size = grid_size
    self.pixels_per_bin = pixels_per_bin

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
    fig = pyplot.figure()
    pyplot.title("Background Model mean histogram")
    pyplot.hist(MEAN, bins=20)
    pyplot.xlabel("mean")
    pyplot.ylabel("# reflections")
    fig.savefig(join(self.directory, "background_model_mean_hist"))
    pyplot.close()

  def mean_vs_xy(self, rlist):
    ''' Plot I/Sigma vs X/Y '''

    class mean_vs_xy_plot(per_panel_plot):

      title = "Distribution of Background Model mean vs X/Y"
      filename = "background_model_mean_vs_xy.png"
      cbar_ylabel = "Background Model mean"

      def plot_one_panel(self, ax, rlist):
        MEAN = rlist['background.mean']
        x, y, z = rlist['xyzcal.px'].parts()

        hex_ax = ax.hexbin(
          x.as_numpy_array(), y.as_numpy_array(),
          C=MEAN.as_numpy_array(), gridsize=self.gridsize,
          vmin=0, vmax=1
        )
        return hex_ax

    plot = mean_vs_xy_plot(rlist, self.directory, grid_size=self.grid_size,
                           pixels_per_bin=self.pixels_per_bin)

  def mean_vs_z(self, rlist):
    ''' Plot I/Sigma vs Z. '''
    from os.path import join
    MEAN = rlist['background.mean']
    x, y, z = rlist['xyzcal.px'].parts()
    fig = pyplot.figure()
    pyplot.title("Distribution of Background Model mean vs Z")
    cax = pyplot.hexbin(z, MEAN, gridsize=100)
    cax.axes.set_xlabel("z")
    cax.axes.set_ylabel("Background Model mean")
    cbar = pyplot.colorbar(cax)
    cbar.ax.set_ylabel("# reflections")
    fig.savefig(join(self.directory, "background_model_mean_vs_z.png"))
    pyplot.close()

  def mean_vs_ios(self, rlist):
    ''' Analyse the correlations. '''
    from os.path import join
    MEAN = rlist['background.mean']
    I = rlist['intensity.sum.value']
    I_sig = flex.sqrt(rlist['intensity.sum.variance'])
    I_over_S = I / I_sig
    mask = I_over_S > 0.1
    I_over_S = I_over_S.select(mask)
    MEAN = MEAN.select(mask)
    fig = pyplot.figure()
    pyplot.title("Background Model mean vs Log I/Sigma")
    cax = pyplot.hexbin(flex.log(I_over_S), MEAN, gridsize=100)
    cbar = pyplot.colorbar(cax)
    cax.axes.set_xlabel("Log I/Sigma")
    cax.axes.set_ylabel("Background Model mean")
    cbar.ax.set_ylabel("# reflections")
    fig.savefig(join(self.directory, "background_model_mean_vs_ios.png"))
    pyplot.close()

  def rmsd_hist(self, rlist):
    ''' Analyse the background RMSD. '''
    from os.path import join
    RMSD = flex.sqrt(rlist['background.mse'])
    MEAN = rlist['background.mean']
    RMSD = RMSD / MEAN
    fig = pyplot.figure()
    pyplot.title("Background Model mean histogram")
    pyplot.hist(RMSD, bins=20)
    pyplot.xlabel("mean")
    pyplot.ylabel("# reflections")
    fig.savefig(join(self.directory, "background_model_cvrmsd_hist"))
    pyplot.close()

  def rmsd_vs_xy(self, rlist):
    ''' Plot I/Sigma vs X/Y '''

    class rmsd_vs_xy_plot(per_panel_plot):

      title = "Distribution of Background Model CVRMSD vs X/Y"
      filename = "background_model_cvrmsd_vs_xy.png"
      cbar_ylabel = "Background Model CVRMSD"

      def plot_one_panel(self, ax, rlist):
        RMSD = flex.sqrt(rlist['background.mse'])
        MEAN = rlist['background.mean']
        RMSD = RMSD / MEAN
        x, y, z = rlist['xyzcal.px'].parts()

        hex_ax = ax.hexbin(
          x.as_numpy_array(), y.as_numpy_array(),
          C=RMSD.as_numpy_array(), gridsize=self.gridsize,
          vmin=0, vmax=1
        )
        return hex_ax

    plot = rmsd_vs_xy_plot(rlist, self.directory, grid_size=self.grid_size,
                           pixels_per_bin=self.pixels_per_bin)

  def rmsd_vs_z(self, rlist):
    ''' Plot I/Sigma vs Z. '''
    from os.path import join
    RMSD = flex.sqrt(rlist['background.mse'])
    MEAN = rlist['background.mean']
    RMSD = RMSD / MEAN
    x, y, z = rlist['xyzcal.px'].parts()
    fig = pyplot.figure()
    pyplot.title("Distribution of Background Model CVRMSD vs Z")
    cax = pyplot.hexbin(z, RMSD, gridsize=100)
    cax.axes.set_xlabel("z")
    cax.axes.set_ylabel("Background Model CVRMSD")
    cbar = pyplot.colorbar(cax)
    cbar.ax.set_ylabel("# reflections")
    fig.savefig(join(self.directory, "background_model_cvrmsd_vs_z.png"))
    pyplot.close()

  def rmsd_vs_ios(self, rlist):
    ''' Analyse the correlations. '''
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
    fig = pyplot.figure()
    pyplot.title("Background Model CVRMSD vs Log I/Sigma")
    cax = pyplot.hexbin(flex.log(I_over_S), RMSD, gridsize=100)
    cbar = pyplot.colorbar(cax)
    cax.axes.set_xlabel("Log I/Sigma")
    cax.axes.set_ylabel("Background Model CVRMSD")
    cbar.ax.set_ylabel("# reflections")
    fig.savefig(join(self.directory, "background_model_cvrmsd_vs_ios.png"))
    pyplot.close()


class IntensityAnalyser(object):
  ''' Analyse the intensities. '''

  def __init__(self, grid_size=None, pixels_per_bin=10):
    from os.path import join

    self.grid_size = grid_size
    self.pixels_per_bin = pixels_per_bin

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
      return {'intensity': {}}

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
      return {'intensity': {}}

    rlist = rlist.select(mask)
    Command.end(" Selected %d integrated reflections" % len(rlist))

    x, y, z = rlist['xyzcal.px'].parts()
    self.nbinsx, self.nbinsy = tuple(
        int(math.ceil(i))
        for i in (flex.max(x)/self.pixels_per_bin,
                  flex.max(y)/self.pixels_per_bin))

    d = OrderedDict()

    # Look at distribution of I/Sigma
    print " Analysing distribution of I/Sigma"
    d.update(self.i_over_s_hist(rlist))
    print " Analysing distribution of I/Sigma vs xy"
    d.update(self.i_over_s_vs_xy(rlist, "sum"))
    print " Analysing distribution of I/Sigma vs xy"
    d.update(self.i_over_s_vs_xy(rlist, "prf"))
    print " Analysing distribution of I/Sigma vs z"
    d.update(self.i_over_s_vs_z(rlist))
    #print " Analysing number of background pixels used"
    #self.num_background_hist(rlist)
    #print " Analysing number of foreground pixels used"
    #self.num_foreground_hist(rlist)

    return {'intensity': d}

  def i_over_s_hist(self, rlist):
    ''' Analyse the correlations. '''
    from os.path import join
    I = rlist['intensity.sum.value']
    I_sig = flex.sqrt(rlist['intensity.sum.variance'])
    I_over_S = I / I_sig
    log_I_over_S = flex.log(I_over_S)
    hist = flex.histogram(log_I_over_S, n_slots=20)

    return {
      'log_i_over_sigma_histogram': {
        'data': [{
          'x': list(hist.slot_centers()),
          'y': list(hist.slots()),
          'type': 'bar',
          'name': 'log_i_over_sigma',
        }],
        'layout': {
          'title': 'Log I/Sigma histogram',
          'xaxis': {'title': 'Log I/Sigma'},
          'yaxis': {'title': 'Number of reflections'},
          'bargap': 0,
        },
      },
    }

  def i_over_s_vs_xy(self, rlist, intensity_type):
    ''' Plot I/Sigma vs X/Y '''

    I_sig = flex.sqrt(rlist['intensity.%s.variance' %intensity_type])
    I = rlist['intensity.%s.value' %intensity_type]
    sel = (I_sig > 0) & (I > 0)
    rlist = rlist.select(sel)
    I = I.select(sel)
    I_sig = I_sig.select(sel)
    I_over_S = I / I_sig
    x, y, z = rlist['xyzcal.px'].parts()

    import numpy as np
    H, xedges, yedges = np.histogram2d(
      x.as_numpy_array(), y.as_numpy_array(), bins=(self.nbinsx, self.nbinsy))
    H1, xedges, yedges = np.histogram2d(
      x.as_numpy_array(), y.as_numpy_array(),
      bins=(self.nbinsx, self.nbinsy),
      weights=flex.log(I_over_S).as_numpy_array())

    nonzeros = np.nonzero(H)
    z = np.empty(H.shape)
    z[:] = np.NAN
    z[nonzeros] = (H1[nonzeros]/H[nonzeros])

    return {
      'i_over_sigma_%s_vs_xy' %intensity_type: {
        'data': [{
          'x': xedges.tolist(),
          'y': xedges.tolist(),
          'z': z.transpose().tolist(),
          'type': 'heatmap',
          'name': 'i_over_sigma_%s' %intensity_type,
          'colorbar': {
            'title': 'Log I/Sigma',
            'titleside': 'right',
          },
          'colorscale': 'Jet',
        }],
        'layout': {
          'title': 'Distribution of I(%s)/Sigma vs X/Y' %intensity_type,
          'xaxis': {
            'domain': [0, 0.85],
            'title': 'X',
            'showgrid': False,
          },
          'yaxis': {
            'title': 'Y',
            'autorange': 'reversed',
            'showgrid': False,
          },
          'width': 500,
          'height': 450,
        },
      },
    }

  def i_over_s_vs_z(self, rlist):
    ''' Plot I/Sigma vs Z. '''

    I = rlist['intensity.sum.value']
    I_sig = flex.sqrt(rlist['intensity.sum.variance'])
    I_over_S = I / I_sig
    x, y, z = rlist['xyzcal.px'].parts()

    import numpy as np
    H, xedges, yedges = np.histogram2d(
      z.as_numpy_array(), flex.log(I_over_S).as_numpy_array(), bins=(100, 100))

    return {
      'i_over_sigma_vs_z': {
        'data': [{
          'x': xedges.tolist(),
          'y': xedges.tolist(),
          'z': H.transpose().tolist(),
          'type': 'heatmap',
          'name': 'i_over_sigma',
          'colorbar': {
            'title': 'Number of reflections',
            'titleside': 'right',
          },
          'colorscale': 'Jet',
        }],
        'layout': {
          'title': 'Difference between observed and calculated centroids in Z',
          'xaxis': {
            'title': 'Z',
            'showgrid': False,
          },
          'yaxis': {
            'title': 'Log I/Sigma',
            'showgrid': False,
          },
        },
      }
    }

  def num_background_hist(self, rlist):
    ''' Analyse the number of background pixels. '''
    from os.path import join
    if 'n_background' in rlist:
      N = rlist['n_background']
      fig = pyplot.figure()
      pyplot.title("Num Background Pixel Histogram")
      pyplot.hist(N, bins=20)
      pyplot.xlabel("Number of pixels")
      pyplot.ylabel("# reflections")
      fig.savefig(join(self.directory, "n_background_hist.png"))
      pyplot.close()

  def num_foreground_hist(self, rlist):
    ''' Analyse the number of foreground pixels. '''
    from os.path import join
    if 'n_foreground' in rlist:
      N = rlist['n_foreground']
      fig = pyplot.figure()
      pyplot.title("Num Foreground Pixel Histogram")
      pyplot.hist(N, bins=20)
      pyplot.xlabel("Number of pixels")
      pyplot.ylabel("# reflections")
      fig.savefig(join(self.directory, "n_foreground_hist.png"))
      pyplot.close()


class ReferenceProfileAnalyser(object):
  ''' Analyse the reference profiles. '''

  def __init__(self, grid_size=None, pixels_per_bin=10):
    from os.path import join

    self.grid_size = grid_size
    self.pixels_per_bin = pixels_per_bin

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
      return {'reference': {}}

    # Select only integrated reflections
    Command.start(" Selecting only integated reflections")
    mask = rlist.get_flags(rlist.flags.integrated)
    if mask.count(True) == 0:
      return {'reference': {}}

    rlist = rlist.select(mask)
    Command.end(" Selected %d integrated reflections" % len(rlist))

    x, y, z = rlist['xyzcal.px'].parts()
    self.nbinsx, self.nbinsy = tuple(
        int(math.ceil(i))
        for i in (flex.max(x)/self.pixels_per_bin,
                  flex.max(y)/self.pixels_per_bin))

    d = OrderedDict()

    # Analyse distribution of reference spots
    print " Analysing reference profile distribution vs x/y"
    d.update(self.reference_xy(rlist))
    print " Analysing reference profile distribution vs z"
    d.update(self.reference_z(rlist))

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

      d = OrderedDict()
      print " Analysing reflection profile correlations"
      d.update(self.reflection_corr_hist(rlist, filename))

      print " Analysing reflection profile correlations vs x/y"
      d.update(self.reflection_corr_vs_xy(rlist, filename))

      print " Analysing reflection profile correlations vs z"
      d.update(self.reflection_corr_vs_z(rlist, filename))

      print " Analysing reflection profile correlations vs I/Sigma"
      d.update(self.reflection_corr_vs_ios(rlist, filename))

      return d

    mask = rlist.get_flags(rlist.flags.reference_spot)
    d.update(correlations("reference",  rlist.select(mask)))
    d.update(correlations("reflection", rlist))
    #ideal_correlations("reference", rlist.select(mask))
    #ideal_correlations("reflection", rlist)

    return {'reference': d}

  def reference_xy(self, rlist):
    ''' Analyse the distribution of reference profiles. '''
    from os.path import join
    mask = rlist.get_flags(rlist.flags.reference_spot)
    rlist = rlist.select(mask)
    x, y, z = rlist['xyzcal.px'].parts()

    import numpy as np
    H, xedges, yedges = np.histogram2d(
      x.as_numpy_array(), y.as_numpy_array(), bins=(self.nbinsx, self.nbinsy))

    return {
      'n_reference_profiles_vs_xy': {
        'data': [{
          'x': xedges.tolist(),
          'y': yedges.tolist(),
          'z': H.transpose().tolist(),
          'type': 'heatmap',
          'name': 'n_reference_profiles',
          'colorbar': {
            'title': 'Number of reflections',
            'titleside': 'right',
          },
          'colorscale': 'Jet',
        }],
        'layout': {
          'title': 'Reference profiles binned in X/Y',
          'xaxis': {
            'domain': [0, 0.85],
            'title': 'X',
            'showgrid': False,
          },
          'yaxis': {
            'title': 'Y',
            'autorange': 'reversed',
            'showgrid': False,
          },
          'width': 500,
          'height': 450,
        },
      },
    }

  def reference_z(self, rlist):
    ''' Analyse the distribution of reference profiles. '''
    corr = rlist['profile.correlation']
    x, y, z = rlist['xyzcal.px'].parts()
    hist = flex.histogram(z, n_slots=20)

    return {
      'n_reference_profiles_vs_z': {
        'data': [{
          'x': list(hist.slot_centers()),
          'y': list(hist.slots()),
          'type': 'bar',
          'name': 'n_reference_profiles',
        }],
        'layout': {
          'title': 'Reference profiles binned in Z',
          'xaxis': {'title': 'Z'},
          'yaxis': {'title': 'Number of reflections'},
          'bargap': 0,
        },
      },
    }

  def reflection_corr_hist(self, rlist, filename):
    ''' Analyse the correlations. '''
    corr = rlist['profile.correlation']
    hist = flex.histogram(corr, n_slots=20)

    return {
      '%s_correlations_histogram' %filename: {
        'data': [{
          'x': list(hist.slot_centers()),
          'y': list(hist.slots()),
          'type': 'bar',
          'name': '%s_correlations' %filename,
        }],
        'layout': {
          'title': '%s correlations histogram' %filename.capitalize(),
          'xaxis': {'title': 'Correlation with reference profile'},
          'yaxis': {'title': 'Number of reflections'},
          'bargap': 0,
        },
      },
    }

  def reflection_corr_vs_xy(self, rlist, filename):
    ''' Analyse the correlations. '''

    corr = rlist['profile.correlation']
    x, y, z = rlist['xyzcal.px'].parts()

    import numpy as np
    H, xedges, yedges = np.histogram2d(
      x.as_numpy_array(), y.as_numpy_array(), bins=(self.nbinsx, self.nbinsy))
    H1, xedges, yedges = np.histogram2d(
      x.as_numpy_array(), y.as_numpy_array(), bins=(self.nbinsx, self.nbinsy),
      weights=corr.as_numpy_array())

    nonzeros = np.nonzero(H)
    z = np.empty(H.shape)
    z[:] = np.NAN
    z[nonzeros] = (H1[nonzeros]/H[nonzeros])

    return {
      '%s_correlations_xy' %filename: {
        'data': [{
          'x': xedges.tolist(),
          'y': yedges.tolist(),
          'z': z.transpose().tolist(),
          'type': 'heatmap',
          'name': '%s_correlations' %filename,
          'colorbar': {
            'title': 'Correlation with reference profile',
            'titleside': 'right',
          },
          'colorscale': 'Jet',
          'zmin': 0,
          'zmax': 1,
        }],
        'layout': {
          'title': '%s correlations binned in X/Y' %filename.capitalize(),
          'xaxis': {
            'domain': [0, 0.85],
            'title': 'X',
            'showgrid': False,
          },
          'yaxis': {
            'title': 'Y',
            'autorange': 'reversed',
            'showgrid': False,
          },
          'width': 500,
          'height': 450,
        },
      },
    }

  def reflection_corr_vs_z(self, rlist, filename):
    ''' Analyse the correlations. '''

    corr = rlist['profile.correlation']
    x, y, z = rlist['xyzcal.px'].parts()

    import numpy as np
    H, xedges, yedges = np.histogram2d(
      z.as_numpy_array(), corr.as_numpy_array(), bins=(100, 100))

    return {
      '%s_correlations_vs_z' %filename: {
        'data': [{
          'x': xedges.tolist(),
          'y': yedges.tolist(),
          'z': H.transpose().tolist(),
          'type': 'heatmap',
          'name': '%s_correlations' %filename,
          'colorbar': {
            'title': 'Number of reflections',
            'titleside': 'right',
          },
          'colorscale': 'Jet',
        }],
        'layout': {
          'title': '%s correlations vs Z' %filename.capitalize(),
          'xaxis': {
            'title': 'Z',
            'showgrid': False,
          },
          'yaxis': {
            'title': 'Correlation with reference profile',
            'showgrid': False,
          },
        },
      },
    }

  def reflection_corr_vs_ios(self, rlist, filename):
    ''' Analyse the correlations. '''

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

    import numpy as np
    H, xedges, yedges = np.histogram2d(
      flex.log(I_over_S).as_numpy_array(),
      corr.as_numpy_array(), bins=(100, 100))

    return {
      '%s_correlations_vs_ios' %filename: {
        'data': [{
          'x': xedges.tolist(),
          'y': yedges.tolist(),
          'z': H.transpose().tolist(),
          'type': 'heatmap',
          'name': '%s_correlations' %filename,
          'colorbar': {
            'title': 'Number of reflections',
            'titleside': 'right',
          },
          'colorscale': 'Jet',
        }],
        'layout': {
          'title': '%s correlations vs Log I/Sigma' %filename.capitalize(),
          'xaxis': {
            'title': 'Log I/Sigma',
            'showgrid': False,
          },
          'yaxis': {
            'title': 'Correlation with reference profile',
            'showgrid': False,
          },
        },
      },
    }

  def ideal_reflection_corr_hist(self, rlist, filename):
    ''' Analyse the correlations. '''
    from os.path import join
    if 'correlation.ideal.profile' in rlist:
      corr = rlist['correlation.ideal.profile']
      pyplot.title("Reflection correlations histogram")
      pyplot.hist(corr, bins=20)
      pyplot.xlabel("Correlation with reference profile")
      pyplot.ylabel("# reflections")
      pyplot.savefig(join(self.directory, "ideal_%s_corr_hist" % filename))
      pyplot.close()

  def ideal_reflection_corr_vs_xy(self, rlist, filename):
    ''' Analyse the correlations. '''
    from os.path import join
    if 'correlation.ideal.profile' in rlist:
      corr = rlist['correlation.ideal.profile']
      x, y, z = rlist['xyzcal.px'].parts()
      pyplot.title("Reflection correlations binned in X/Y")
      cax = pyplot.hexbin(x, y, C=corr, gridsize=100, vmin=0.0, vmax=1.0)
      cbar = pyplot.colorbar(cax)
      pyplot.xlabel("x")
      pyplot.ylabel("y")
      cbar.ax.set_ylabel("Correlation with reference profile")
      pyplot.savefig(join(self.directory, "ideal_%s_corr_vs_xy.png" % filename))
      pyplot.close()

  def ideal_reflection_corr_vs_z(self, rlist, filename):
    ''' Analyse the correlations. '''
    from os.path import join
    if 'correlation.ideal.profile' in rlist:
      corr = rlist['correlation.ideal.profile']
      x, y, z = rlist['xyzcal.px'].parts()
      pyplot.title("Reflection correlations vs Z")
      cax = pyplot.hexbin(z, corr, gridsize=100)
      cbar = pyplot.colorbar(cax)
      pyplot.xlabel("z")
      pyplot.ylabel("Correlation with reference profile")
      cbar.ax.set_ylabel("# reflections")
      pyplot.savefig(join(self.directory, "ideal_%s_corr_vs_z.png" % filename))
      pyplot.close()

  def ideal_reflection_corr_vs_ios(self, rlist, filename):
    ''' Analyse the correlations. '''
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
      pyplot.title("Reflection correlations vs Log I/Sigma")
      cax = pyplot.hexbin(flex.log(I_over_S), corr, gridsize=100)
      cbar = pyplot.colorbar(cax)
      pyplot.xlabel("Log I/Sigma")
      pyplot.ylabel("Correlation with reference profile")
      cbar.ax.set_ylabel("# reflections")
      pyplot.savefig(join(self.directory, "ideal_%s_corr_vs_ios.png" % filename))
      pyplot.close()


class Analyser(object):
  ''' Helper class to do all the analysis. '''

  def __init__(self, params, grid_size=None, centroid_diff_max=1.5):
    ''' Setup the analysers. '''
    self.params = params
    self.analysers = [
      StrongSpotsAnalyser(),
      CentroidAnalyser(
        grid_size=grid_size, pixels_per_bin=self.params.pixels_per_bin,
          centroid_diff_max=centroid_diff_max),
      BackgroundAnalyser(
         grid_size=grid_size, pixels_per_bin=self.params.pixels_per_bin),
      IntensityAnalyser(
        grid_size=grid_size, pixels_per_bin=self.params.pixels_per_bin),
      ReferenceProfileAnalyser(
        grid_size=grid_size, pixels_per_bin=self.params.pixels_per_bin),
    ]

  def __call__(self, rlist=None, experiments=None):
    ''' Do all the analysis. '''
    from copy import deepcopy
    json_data = OrderedDict()

    if rlist is not None:
      for analyse in self.analysers:
        result = analyse(deepcopy(rlist))
        if result is not None:
          json_data.update(result)
    else:
      json_data.update(
        {'strong': {}, 'centroid': {}, 'intensity': {}, 'reference': {}})

    if experiments is not None:
      analyse = ScanVaryingCrystalAnalyser(self.params.orientation_decomposition)
      json_data.update(analyse(experiments))

    import json
    json_str = json.dumps(json_data)


    if self.params.output.html is not None:
      javascript = ['var graphs = %s' %(json_str)]

      graph_divs = {}
      for grouping in json_data.keys():
        graph_divs[grouping] = []
        for graph in json_data[grouping].keys():
          javascript.append(
            'Plotly.newPlot(%(graph)s, graphs.%(grouping)s.%(graph)s.data, graphs.%(grouping)s.%(graph)s.layout);' %{
              'graph': graph,
              'grouping': grouping
            })

          graph_divs[grouping].append(
            '<div class="col-xs-6 col-sm-6 col-md-4 plot" id="%(graph)s"></div>' %{'graph': graph})

      html_header = '''
<head>

<!-- Plotly.js -->
<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>

<meta name="viewport" content="width=device-width, initial-scale=1" charset="UTF-8">
<link rel="stylesheet" href="http://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/css/bootstrap.min.css">
<script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.3/jquery.min.js"></script>
<script src="http://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/js/bootstrap.min.js"></script>

<style type="text/css">

body {
  /*font-family: Helmet, Freesans, Helvetica, Arial, sans-serif;*/
  margin: 8px;
  min-width: 240px;
  margin-left: 5%;
  margin-right: 5%;
}

.plot {
  float: left;
  width: 600px;
  height: 400px;
  margin-bottom: 20px;
}

</style>

</head>

'''

      html_body = '''

<body>

<div class="page-header">
  <h1>DIALS analysis report</h1>
</div>

<div >
  <h2>DIALS analysis plots</h2>
  <div class="panel-group">
    <div class="panel panel-default">
      <div class="panel-heading" data-toggle="collapse" href="#collapse1">
        <h4 class="panel-title">
          <a>Analysis of strong reflections</a>
        </h4>
      </div>
      <div id="collapse1" class="panel-collapse collapse">
        <div class="panel-body">

          <div class="container-fluid">
            %(strong_graph_divs)s
          </div>

        </div>
        <!-- <div class="panel-footer"></div> -->
      </div>
    </div>
    <div class="panel panel-default">
      <div class="panel-heading" data-toggle="collapse" href="#collapse2">
        <h4 class="panel-title">
          <a>Analysis of reflection centroids</a>
        </h4>
      </div>
      <div id="collapse2" class="panel-collapse collapse">
        <div class="panel-body">

          <div class="container-fluid">
            %(centroid_graph_divs)s
          </div>

        </div>
        <!-- <div class="panel-footer"></div> -->
      </div>
    </div>
    <div class="panel panel-default">
      <div class="panel-heading" data-toggle="collapse" href="#collapse3">
        <h4 class="panel-title">
          <a>Analysis of reflection intensities</a>
        </h4>
      </div>
      <div id="collapse3" class="panel-collapse collapse">
        <div class="panel-body">

          <div class="container-fluid">
            %(intensity_graph_divs)s
          </div>

        </div>
        <!-- <div class="panel-footer"></div> -->
      </div>
    </div>
    <div class="panel panel-default">
      <div class="panel-heading" data-toggle="collapse" href="#collapse4">
        <h4 class="panel-title">
          <a>Analysis of reference profiles</a>
        </h4>
      </div>
      <div id="collapse4" class="panel-collapse collapse">
        <div class="panel-body">

          <div class="container-fluid">
            %(reference_graph_divs)s
          </div>

        </div>
        <!-- <div class="panel-footer"></div> -->
      </div>
    </div>
    <div class="panel panel-default">
      <div class="panel-heading" data-toggle="collapse" href="#collapse5">
        <h4 class="panel-title">
          <a>Analysis of scan-varying crystal model</a>
        </h4>
      </div>
      <div id="collapse5" class="panel-collapse collapse">
        <div class="panel-body">

          <div class="container-fluid">
            %(scan_varying_graph_divs)s
          </div>

        </div>
        <!-- <div class="panel-footer"></div> -->
      </div>
    </div>
  </div>
</div>

<script>
%(script)s
</script>
</body>
''' %{'strong_graph_divs': '\n            '.join(graph_divs['strong']),
      'centroid_graph_divs': '\n            '.join(graph_divs['centroid']),
      'intensity_graph_divs': '\n            '.join(graph_divs['intensity']),
      'reference_graph_divs': '\n            '.join(graph_divs['reference']),
      'scan_varying_graph_divs': '\n            '.join(graph_divs['scan_varying']),
      'script': '\n'.join(javascript)}

      html = '\n'.join([html_header, html_body])

      print "Writing html report to: %s" %self.params.output.html
      with open(self.params.output.html, 'wb') as f:
        print >> f, html.encode('ascii', 'xmlcharrefreplace')

    if self.params.output.json is not None:
      print "Writing json data to: %s" %self.params.output.json
      with open(self.params.output.json, 'wb') as f:
        print >> f, json_str


class Script(object):
  ''' A class to encapsulate the script. '''

  def __init__(self):
    ''' Initialise the script. '''
    from dials.util.options import OptionParser
    import libtbx.load_env

    # Create the parser
    usage = "usage: %s [options] reflections.pickle" % libtbx.env.dispatcher_name
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      read_reflections=True,
      read_experiments=True,
      epilog=help_message)

  def run(self):
    ''' Run the script. '''
    from dials.util.command_line import Command
    from libtbx.utils import Sorry
    from dials.util.options import flatten_reflections, flatten_experiments

    # Parse the command line arguments
    params, options = self.parser.parse_args(show_diff_phil=True)

    # Shoe the help
    if len(params.input.reflections) != 1 and not len(params.input.experiments):
      self.parser.print_help()
      exit(0)

    from dials.util.options import flatten_reflections, flatten_experiments
    reflections = flatten_reflections(params.input.reflections)
    experiments = flatten_experiments(params.input.experiments)

    # Analyse the reflections
    analyse = Analyser(params, grid_size=params.grid_size,
      centroid_diff_max=params.centroid_diff_max)
    if len(reflections):
      reflections = reflections[0]
    else:
      reflections = None

    analyse(reflections, experiments)


if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
