#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2014) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

# LIBTBX_SET_DISPATCHER_NAME dev.dials.plot_detector_shifts

"""
Given two detector models (e.g. one refined at hierarchy_level=0 and another
at hierarchy_level=1), plot shifts in the fast, slow and normal directions on
each panel that take the first detector model to the second one as a heatmap

"""
from __future__ import division

from dxtbx.model.experiment.experiment_list import Experiment
from scitbx.array_family import flex
from libtbx import phil, Auto
from libtbx.utils import Sorry
import matplotlib
# Offline backend
matplotlib.use("Agg")
matplotlib.rc('font', family='serif')
matplotlib.rc('font', serif='Times New Roman')
from matplotlib import pyplot as plt

help_message = '''

Plots shifts between two detector models. The shifts are those that take pixel
positions from the first detector (the reference) to the equivalent position
on the second detector, and are expressed in the basis (fast, slow, normal) of
each panel of the reference detector.

Example::

  dev.dials.plot_detector_shifts level0.json level1.json

Here level0 might contain a hierarchical detector refined at hierarchy_level=0
(i.e. bulk movements), while level1.json may contain a detector with metrology
refined one level below that, starting from level0.json. The resulting plot
will show the shifts that occurred during the refinement that resulted in
level1.json.
'''

phil_scope = phil.parse('''
plot_type = *panel_grid spherical_polar
  .type = choice
  .help = "Choose type of plot"
grid_size = Auto
  .type = ints(size=2)
  .help = "Grid size (nrow, ncol) for panel_grid plots. Set automatically for"
          "known detector types if set to Auto"
tag = None
  .type = str
  .help = "output files will be pre-pended with this string"
''')

# sample 1 pt per mm
SAMPLE_FREQ = 1

class PlotData(object):

  def __init__(self, detector1, detector2):

    self.det1 = detector1
    self.det2 = detector2

  def __call__(self, ipanel=0):

    panel_a = self.det1[ipanel]
    panel_b = self.det2[ipanel]
    size_fast, size_slow = panel_a.get_image_size_mm()
    assert size_fast, size_slow == panel_b.get_image_size_mm()

    # num of sample intervals
    n_fast = int((size_fast) / SAMPLE_FREQ)
    n_slow = int((size_slow) / SAMPLE_FREQ)

    # interval width
    step_fast = size_fast / n_fast
    step_slow = size_slow / n_slow

    # samples
    samp_fast = [step_fast * i for i in range(n_fast + 1)]
    samp_slow = [step_slow * i for i in range(n_slow + 1)]

    lab1 = flex.vec3_double()
    lab2 = flex.vec3_double()
    sample_pts = flex.vec2_double()

    # loop
    for s in samp_slow:
      for f in samp_fast:
        lab1.append(panel_a.get_lab_coord((f, s)))
        lab2.append(panel_b.get_lab_coord((f, s)))
        sample_pts.append((f,s))

    offset = lab2 - lab1

    # store offset in the lab frame
    x_off, y_off, z_off = offset.parts()

    # reexpress offset in the basis fast, slow, normal of panel_a
    f_off = offset.dot(panel_a.get_fast_axis())
    s_off = offset.dot(panel_a.get_slow_axis())
    n_off = offset.dot(panel_a.get_normal())

    f, s = sample_pts.parts()

    return {'lab_coord':lab1,
            'fast':f, 'slow':s,
            'x_offset':x_off,
            'y_offset':y_off,
            'z_offset':z_off,
            'fast_offset':f_off,
            'slow_offset':s_off,
            'normal_offset':n_off,
            'size_fast':size_fast,
            'size_slow':size_slow}

def plot_grid_of_panels(panel_data, nrow, ncol, direction='fast', tag = ''):
  '''Plot data for each panel in a stack of subplots, with the first panel
  at the top. This is appropriate for e.g. the 24 panel model for the I23
  P12M detector, as it is simply unrolling the barrel to make a flat plot'''

  fig, axarr = plt.subplots(nrow, ncol, sharex=True)
  # kludge just in case nrow=ncol=1 and axarr is not actually an array
  try:
    len(axarr)
  except TypeError:
    import numpy as np
    axarr = np.array(axarr)
  plt.suptitle(r"$\Delta " + direction + "$" + " shifts", y=0.95)
  plt.figtext(0.425, 0.05, "fast (mm)")

  plt.setp(axarr.flat, aspect=1.0, adjustable='box-forced')

  clim = []
  imarr = []
  for ipanel, (pnl_data, ax) in enumerate(zip(panel_data, axarr.flatten())):
    f, s, offset = (pnl_data['fast'], pnl_data['slow'],
                    pnl_data[direction + '_offset'])
    im = ax.hexbin(f.as_numpy_array(),
                           s.as_numpy_array(),
                           offset.as_numpy_array(),
                           gridsize=30)
    imarr.append(im)
    clim.append(im.get_clim())
    ax.invert_yaxis()
    ax.set_yticks([]) # don't show axis labels on the stacked side
    ax.set_xticks([0, round(pnl_data['size_fast'])])
    ax.tick_params('y', labelsize='x-small')

  clim_low, clim_high = zip(*clim)
  clim = (min(clim_low), max(clim_high))
  if (clim[1] - clim[0]) < 1e-12:
    print "...skipping plot with shift too small to show"
    return
  for im in imarr: im.set_clim(clim)
  fig.subplots_adjust(right=0.8)
  cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
  fig.colorbar(im, cax=cbar_ax)
  cbar_ax.set_ylabel("$" + direction + "_{2} - " + direction + "_{1}$ (mm)")
  fig.set_size_inches(10,10)
  plt.savefig(tag + direction + "_diff.png")
  plt.clf()

def plot_spherical_polar(panel_data, beam, direction='fast', tag = ''):
  '''Plot data for all panels in a single plot by mapping pixel positions to
  the surface of the Ewald sphere, then plotting in 2D using azimuth and
  elevation angles. This distorts the image from a flat panel detector, but
  will work for any detector geometry. The equatorial plane for the spherical
  coordinate system is defined by a beam direction and the lab X direction.
  The azimuth and elevation angles are zero along the beam direction.'''

  # lists for the plot axes
  azimuth = flex.double()
  elevation = flex.double()
  offset = flex.double()

  # define equatorial plane with two vectors, s0u and locx
  from scitbx import matrix
  s0u = matrix.col(beam.get_unit_s0())
  labx = matrix.col((1, 0, 0))
  norm = labx.cross(s0u).normalize()
  locx = s0u.cross(norm).normalize()

  # calculate components of the lab coords in the equatorial plane directions
  for dat in panel_data:
    lab = dat['lab_coord']
    off = dat[direction + '_offset']
    a = lab.dot(s0u)
    b = lab.dot(locx)

    # calculate vectors in the equatorial plane
    plane_proj = a * flex.vec3_double(len(a), s0u) + \
                 b * flex.vec3_double(len(b), locx)

    # hence azimuthal angles
    az = plane_proj.angle(s0u, deg=True)
    neg = plane_proj.dot(locx) < 0.0
    az.set_selected(neg, -1.0 * az.select(neg))

    # and elevation angles
    el = lab.angle(plane_proj, deg=True)
    neg = lab.dot(norm) < 0.0
    el.set_selected(neg, -1.0 * el.select(neg))

    azimuth.extend(az)
    elevation.extend(el)
    offset.extend(off)

  fig=plt.figure()
  plt.xlabel("azimuth (degrees)")
  plt.ylabel("elevation (degrees)")
  ax=fig.add_subplot(111)
  ax.set_title(r"$\Delta " + direction + "$")
  im = ax.hexbin(azimuth.as_numpy_array(),
                 elevation.as_numpy_array(),
                 offset.as_numpy_array(),
                 gridsize=100)
  ax.set_xlim(flex.min(azimuth)*1.1, flex.max(azimuth)*1.1)
  ax.set_ylim(flex.min(elevation)*1.1, flex.max(elevation)*1.1)
  fig.subplots_adjust(right=0.8)
  cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
  fig.colorbar(im, cax=cbar_ax)
  cbar_ax.set_ylabel("$" + direction + "_{2} - " + direction + "_{1}$ (mm)")
  plt.savefig(tag + direction + "_diff.png")
  plt.clf()

  return

def guess_grid_size(detector):
  '''Guess the right grid size for the detector based on comparison with
  known types'''

  n_panels = len(detector)
  pxl_sizes = [p.get_pixel_size() for p in detector]
  pnl_sizes = [p.get_image_size() for p in detector]

  # check all sizes are equal
  if not all([s == pxl_sizes[0] for s in pxl_sizes[1:]]):
    return None
  if not all([s == pnl_sizes[0] for s in pnl_sizes[1:]]):
    return None

  # take the first
  pxl_size = [round(s, 3) for s in pxl_sizes[0]]
  pnl_size = list(pnl_sizes[0])

  if n_panels == 24 and pxl_size == [0.172, 0.172] and pnl_size == [2463, 195]:
    # guessed P12M 24 row model
    print "Guessed P12M detector modelled with 24 rows"
    grid_size = (24, 1)
  elif n_panels == 120 and pxl_size == [0.172, 0.172] and pnl_size == [487, 195]:
    # guessed P12M 24 row, 5 column model
    print "Guessed P12M detector modelled with 24 rows in 5 columns"
    grid_size = (24, 5)
  elif n_panels == 60 and pxl_size == [0.172, 0.172] and pnl_size == [487, 195]:
    print "Guessed P6M detector modelled with 60 separate panels"
    grid_size = (12, 5)
  elif n_panels == 64 and pxl_size == [0.110, 0.110] and pnl_size == [194, 185]:
    print "Guessed CS-PAD detector. WARNING plot_type=spherical_polar is probably more useful!"
    grid_size = (8, 8)
  else:
    grid_size = None

  return grid_size

class Script(object):

  def __init__(self):
    '''Check script input and return two experiments if all is okay'''

    import libtbx.load_env
    from dials.util.options import OptionParser

    # The script usage
    usage = ("usage: {0} [options] [param.phil] experiments1.json "
             "experiments2.json").format(libtbx.env.dispatcher_name)

    parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      read_datablocks=True,
      read_experiments=True,
      check_format=False,
      epilog=help_message)

    params, options = parser.parse_args(show_diff_phil=True)

    if params.input.datablock:
      if params.input.experiments:
        raise Sorry("Please provide either datablocks or experiment lists, "
                    "not a combination of the two")
      if len(params.input.datablock) != 2:
        raise Sorry("Please provide two datablocks or experiment lists as input")

      warnmsg = ("WARNING: The {0} datablock contains more than one "
                 "detector. Only the first will be considered.")
      detector1 = params.input.datablock[0].data[0].unique_detectors()
      if len(detector1) > 1:
        print warnmsg.format("first")
      detector1 = detector1[0]
      beam1 = params.input.datablock[0].data[0].unique_beams()[0]
      experiment1 = Experiment(beam=beam1, detector=detector1)

      detector2 = params.input.datablock[1].data[0].unique_detectors()
      if len(detector2) > 1:
        print warnmsg.format("second")
      detector2 = detector2[0]
      beam2 = params.input.datablock[1].data[0].unique_beams()[0]
      experiment2 = Experiment(beam=beam2, detector=detector2)

    else:
      if len(params.input.experiments) != 2:
        raise Sorry("Please provide two datablocks or experiment lists as input")

      warnmsg = ("WARNING: The {0} experiment list contains more than one "
                 "detector. Only the first will be considered.")
      detector1 = params.input.experiments[0].data.detectors()
      if len(detector1) > 1:
        print warnmsg.format("first")
      detector1 = detector1[0]
      experiment1 = params.input.experiments[0].data[0]

      detector2 = params.input.experiments[1].data.detectors()
      if len(detector2) > 1:
        print warnmsg.format("second")
      detector2 = detector2[0]
      experiment2 = params.input.experiments[1].data[0]

    if len(detector1) != len(detector2):
      raise Sorry("The detectors do not contain the same number of panels")

    self.experiment1 = experiment1
    self.experiment2 = experiment2
    self.params = params
    return

  def __call__(self):

    det1, det2 = self.experiment1.detector, self.experiment2.detector

    plot_data = PlotData(det1, det2)

    # do calculations in advance and store in a dictionary
    dat=[]
    for ipanel in range(len(det1)):
      print "Calc for panel", ipanel
      dat.append(plot_data(ipanel))

    if self.params.tag is None:
      tag = ''
    else:
      tag = "%s_"%self.params.tag

    grid_size = self.params.grid_size
    if self.params.plot_type == 'panel_grid' and grid_size is Auto:
      grid_size = guess_grid_size(det1)
      if grid_size is None:
        raise Sorry("Unable to guess the grid_size for this detector")

    # first the plots using offsets in local fast, slow, normal frames
    for direction in ['fast', 'slow', 'normal']:
      print "Doing plot of offsets in the {0} direction".format(direction)
      if self.params.plot_type == 'panel_grid':
        plot_grid_of_panels(dat, grid_size[0], grid_size[1], direction, tag)
      elif self.params.plot_type == 'spherical_polar':
        plot_spherical_polar(dat, self.experiment1.beam, direction, tag)

    # now plots using offsets in the laboratory frame
    for direction in ['x', 'y', 'z']:
      print "Doing plot of offsets in the laboratory {0} direction".format(direction)
      if self.params.plot_type == 'panel_grid':
        plot_grid_of_panels(dat, grid_size[0], grid_size[1], direction, tag)
      elif self.params.plot_type == 'spherical_polar':
        plot_spherical_polar(dat, self.experiment1.beam, direction, tag)

if __name__ == "__main__":

  run = Script()
  run()
