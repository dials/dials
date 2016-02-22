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

import sys
from dxtbx.model.experiment.experiment_list import ExperimentListFactory
import matplotlib
# Offline backend
matplotlib.use("Agg")
matplotlib.rc('font', family='serif')
matplotlib.rc('font', serif='Times New Roman')
from matplotlib import pyplot as plt
from scitbx.array_family import flex

# sample 1 pt per mm
SAMPLE_FREQ = 1

def get_plot_data(ipanel=0):

  # panel 0 first
  panel_a = det1[ipanel]
  panel_b = det2[ipanel]
  size_fast, size_slow = panel_a.get_image_size_mm()
  assert size_fast, size_slow == panel_b.get_image_size_mm()
  aspect = max(size_fast, size_slow) / min(size_fast, size_slow)

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

  # reexpress offset in the basis fast, slow, normal of panel_a
  f_off = offset.dot(panel_a.get_fast_axis())
  s_off = offset.dot(panel_a.get_slow_axis())
  n_off = offset.dot(panel_a.get_normal())

  f, s = sample_pts.parts()

  return f, s, f_off, s_off, n_off


if __name__ == "__main__":

  assert len(sys.argv) == 3

  exp1 = ExperimentListFactory.from_json_file(sys.argv[1],
                check_format=False)
  exp2 = ExperimentListFactory.from_json_file(sys.argv[2],
                check_format=False)

  # take the first experiment only
  det1 = exp1[0].detector
  det2 = exp2[0].detector

  assert len(det1) == len(det2)

  # do calculations in advance and store in a dictionary
  dat={}
  for ipanel in range(len(det1)):
    print "Calc for panel", ipanel
    f, s, f_off, s_off, n_off = get_plot_data(ipanel)
    dat[ipanel] = f, s, f_off, s_off, n_off

  # set limits on colour scale to shifts of 2 pixels
  extrema = 1*0.172

  # first the fast plot
  print "Doing plot of offsets in the fast direction"
  fig, axarr = plt.subplots(len(det1), sharex=True)
  axarr[0].set_title(r"$\Delta fast$")
  plt.xlabel("fast (mm)")

  for ipanel in range(len(det1)):
    f, s, f_off, s_off, n_off = dat[ipanel]
    ax=axarr[ipanel]
    im = ax.hexbin(list(f), list(s), C=list(f_off), gridsize=30)
    ax.invert_yaxis()
    ax.set_yticks([0., 20.])
    ax.tick_params('y', labelsize='x-small')

  fig.subplots_adjust(right=0.8)
  cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
  fig.colorbar(im, cax=cbar_ax)
  cbar_ax.set_ylabel(r"$fast_{2} - fast_{1}$ (mm)")
  plt.savefig("fast_diff.png")
  plt.clf()

  # now the slow plot
  print "Doing plot of offsets in the slow direction"
  fig, axarr = plt.subplots(len(det1), sharex=True)
  axarr[0].set_title(r"$\Delta slow$")
  plt.xlabel("fast (mm)")

  for ipanel in range(len(det1)):
    f, s, f_off, s_off, n_off = dat[ipanel]
    ax=axarr[ipanel]
    im = ax.hexbin(list(f), list(s), C=list(s_off), gridsize=30)
    ax.invert_yaxis()
    ax.set_yticks([0., 20.])
    ax.tick_params('y', labelsize='x-small')

  fig.subplots_adjust(right=0.8)
  cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
  fig.colorbar(im, cax=cbar_ax)
  cbar_ax.set_ylabel(r"$slow_{2} - slow_{1}$ (mm)")
  plt.savefig("slow_diff.png")
  plt.clf()

  # finally the normal plot
  print "Doing plot of offsets in the normal direction"
  fig, axarr = plt.subplots(len(det1), sharex=True)
  axarr[0].set_title(r"$\Delta normal$")
  plt.xlabel("fast (mm)")

  for ipanel in range(len(det1)):
    f, s, f_off, s_off, n_off = dat[ipanel]
    ax=axarr[ipanel]
    im = ax.hexbin(list(f), list(s), C=list(n_off), gridsize=30)
    ax.invert_yaxis()
    ax.set_yticks([0., 20.])
    ax.tick_params('y', labelsize='x-small')

  fig.subplots_adjust(right=0.8)
  cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
  fig.colorbar(im, cax=cbar_ax)
  cbar_ax.set_ylabel(r"$normal_{2} - normal_{1}$ (mm)")
  plt.savefig("normal_diff.png")
  plt.clf()
